####################################################
#### Filtering the Petrov's Diffusion - Library ####
####################################################

# This file contains the low-level implementation of the methods used in the
# manuscript.
#
# In this directory:
# - documented project setup lives in `R/pdhmm_setup.R`;
# - a reusable, documented public interface lives in `R/pdhmm_api.R`;
# - manuscript-specific workflows live in `analysis/`.
#
# The code below is organised in the following broad blocks:
# 1. numerical helpers and scoring rules;
# 2. Poisson--Dirichlet sampling and posterior simulation;
# 3. ancestral-process and diffusion simulation utilities;
# 4. partition algebra and coagulation operators;
# 5. update and prediction recursions;
# 6. particle-filter and predictive wrappers;
# 7. bootstrap particle filter;
# 8. smoothing and cost-to-go routines.

library(DescTools) #Combinations and Permutations
library(partitions)
library(collections) #Dictionary to store coagulations (keys) and their probabilities (values)
library(LaplacesDemon) #Dirichlet distribution
library(extraDistr) #multivariate hypergeometric
library(scales)
library(vioplot)
library(pbapply)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(RcppAlgos)
library(microbenchmark)
library(dplyr)
library(pbmcapply)
library(R.utils)
library(geomtextpath)
library(igraph)
library(lubridate)

#--------------------------------------------------

FallingFactorial <- function(a, b, log = FALSE) {
  if (a < b) return(0)
  if (a == 0) {
    if(log) return(0) else return(1)
  }
  if (log) return(sum(log(a - 0:(b-1) ))) else return( prod( a-0:(b-1) ) )
}

RisingFactorial <- function(a, b, log = FALSE) {
  if (b == 0) {
    if (log) return(0) else return(1)
  }
  
  if (a+b < 0) {
    if (log) return(NULL) else return(gamma(a + b) / gamma(a)) #can't do log if a + b < 0 because gamma(a+b) < 0
  }

  if (a+b>=0 && b < 0) {
    if (log) return(log(gamma(a+b)) - log(gamma(a)) ) else return(gamma(a + b) / gamma(a)) 
  }
  
  if (log) return(sum(log(a + 0:(b-1) ))) else return( prod( a+0:(b-1) ) )
}

#--------------------------------------------------

qentropy <- function(p, q) {
  if (q == 0) return(length(p))
  
  if (q == 1) {
    p <- p[p > 0]
    return( -sum( p * log(p) ) ) 
  }
  else {
    return( (1 - sum(p^q)) / (q - 1) )
  }
}

#--------------------------------------------------
#Interval score of Gneiting and Raftery 

interval_score <- function(qs, p, alpha = 0.05) {
  scores <- c()
  for (i in 1:length(p)) {
    temp <- qs[2, i] - qs[1, i]
    if (p[i] < qs[1, i]) temp <- temp + 2*(qs[1, i] - p[i]) / alpha
    if (p[i] > qs[2, i]) temp <- temp + 2*(p[i] - qs[2, i]) / alpha
    scores <- c(scores, temp)
  }
  return( mean( scores ) )
}

coverage <- function(qs, p) {
  temp1 <- p > qs[1,]
  temp2 <- p < qs[2,]
  temp <- temp1 * temp2
  return(sum(temp) / length(p))
}

#--------------------------------------------------
# epsilon-truncated Poisson-Dirichlet sampling

gem <- function(alpha, theta, epsilon, min_length = 0){
  w <- c()
  residual <- 1
  stick <- 1
  
  i <- 0
  
  while(residual > epsilon || i < min_length){
    i <- i + 1
    
    v_i <- rbeta(1, 1 - alpha, theta + i*alpha)
    w_i <- v_i * stick
    
    stick <- stick * (1-v_i)
    w <- c(w, w_i)
    residual <- residual - w_i
  }
  
  return(list("weights" = w, "residual" = residual))
}

pd <- function(alpha, theta, epsilon, show_top = Inf, min_length = 0){
  w <- gem(alpha, theta, epsilon, min_length = min_length)
  
  x <- sort(w$weights, decreasing = TRUE)
  residual <- w$residual
  
  if ( show_top < Inf && length(x) < show_top ) x <- c(x, rep(0, show_top - length(x) ))
  if (show_top < Inf && length(x) > show_top) x <- x[1:show_top]
  
  return( list("weights" = x, "residual" = residual) )
}

pd_posterior <- function(alpha, theta, partition, epsilon, show_top = Inf, min_length = 0){
  
  m <- sum(partition) #; print(paste('m', m))
  r <- length(partition) #; print(paste('r', length(partition)))
  #print(paste('partition', partition))
  
  #print(paste('BETA', m - r*alpha, theta + r*alpha ))
  if(m == 0) Z <- 0 else Z <- rbeta(1, m - r*alpha, theta + r*alpha)
  #print(Z)
  
  if(m == 0) dir_component <- 0
  if(r == 1) dir_component <- 1 else dir_component <- rdirichlet(1, partition - alpha)
  #degenerate dirichlet might happen if length(partition) = 1
  
  #print(paste('Z', Z, epsilon / (1 - Z) ) )
  pd_component <- pd(alpha, theta + r*alpha, epsilon = epsilon / (1 - Z), min_length = max(0, min_length - r ) )$weights
  
  mixture <- c(Z*dir_component, (1-Z)*pd_component)
  mixture <- sort(mixture, decreasing = T)
  
  if ( show_top < Inf && length(mixture) < show_top ) mixture <- c(mixture, rep(0, show_top - length(mixture) ))
  if (show_top < Inf && length(mixture) > show_top) mixture <- mixture[1:show_top]
  
  return(list('weights'=mixture, 'Dir'=dir_component, 'PD' = pd_component) )
  
}

#Ewens-Pitman Sampling Formula
PSF <- function( alpha, theta, partition, normalized = FALSE ) { #pitman sampling formula
  if ( prod(partition == c(0)) == 1 ) return(1)
  
  L <- length(partition)
  
  cum_partition <- c(0, cumsum(partition))
  
  factorization <- c()
  
  for (l in 1:L){
    V <- (l-1) * alpha + theta
    V <- if (V > 0) V else 1
    if (partition[l] > 1) W <- (1 - alpha ) + 0:(partition[l]  - 2) else W <- 1
    if (partition[l] > 1) temp_num <- c(V, W) else temp_num <- V
    #print(temp_num)
    temp_den <- theta + cum_partition[l]:(cum_partition[l+1] - 1)
    temp_den[temp_den == 0] <- 1
    #print(temp_den)
    temp <- temp_num / temp_den
    #print(temp)
    factorization <- c(factorization, prod(temp))
  }
  
  #print(factorization)
  prob <- prod(factorization)
  
  if( is.nan(prob) ) print(paste('WARNING! NaN produced in PSF! Partition = ', toString(partition), ', parameters = ', toString(c(alpha, theta) ) ) )
  
  if(normalized) return( C(partition)*prob ) else return(prob)
  
}

#--------------------------------------------------
# Simulating the 2-par. Poisson-Dirichlet diffusion (a.k.a. Petrov's diffusion)

#Kingman's paintbox
paintbox <- function(n, x){
  
  if (n == 0) return(c(0))
  
  x <- x / sum(x)
  cumulative <- cumsum(x)
  S <- sum(x)
  
  U <- runif(n)
  
  #count_new <- 0 #use this if yoou want the missing mass to be 'dust'
  
  #pi <- sapply(U, function(u) if(u <= S) which.max(u <= cumulative) else count_new <- count_new + 1) #use this if yoou want the missing mass to be 'dust'
  pi <- sapply(U, function(u) which.max(u <= cumulative) )
  
  #partition <- c(as.numeric(sort(table(pi), decreasing = T)), rep(1, count_new) ) #use this if yoou want the missing mass to be 'dust'
  partition <- as.numeric(sort(table(pi), decreasing = T))
  
  return(partition)
  
}
  
#--------------------------------------------------
# Exact simulation of the ancestral process (Jenkins, Spanò 2017)

a_km <- function(k, m, theta){
  if (k < m) return(NULL)
  
  if (m == 0 && k == 0) return(1)
  #by continuity: (theta - 1) / (theta - 1) -> 1, although in theta = 0 is not defined
  #if return(0) when theta = 0, Ancestral never converges...
  
  num1 <- theta + 2 * k - 1
  
  if (num1 < 0) {
    num2 <- RisingFactorial(theta + m, k - 1, log = F)
    den1 <- FallingFactorial(m, m, log = T)
    den2 <- FallingFactorial(k-m, k-m, log = T)
    den <- den1 + den2
    
    return(num1*num2*exp(-den))
  }
  
  num2 <- RisingFactorial(theta + m, k - 1, log = T)
  #problematic when theta = 1 and m = k = 0 since gamma(0) is Inf!
  #solved by returning 0 since num1 = 0...better ask Jenkins, Spano if it is right...
  den1 <- FallingFactorial(m, m, log = T)
  den2 <- FallingFactorial(k-m, k-m, log = T)
  den <- den1 + den2
  
  return(num1*exp(num2 - den))
}

b <- function(k, m, t, theta){
  
  fact1 <- a_km(k, m, theta)
  
  fact2 <- exp(- ( k * (k + theta - 1) * t ) / 2 )
  
  if( is.na(fact1) ) {print('WARNING: NaN produced in a_km!') ; print(paste(k, m, t, theta)) }
  if( is.na(fact2) ) print('WARNING: NaN produced in b!')
  
  return( fact1 * fact2 )
  
}

b_frac <- function(k, m, t, theta) {
  fact1 <- (theta + 2*k + 1) / (theta + 2*k - 1)
  fact2 <- (theta + m + k - 1) / (k+1-m)
  fact3 <- exp( - (2*k + theta) * t / 2)
  
  return(fact1 * fact2 * fact3)
}
  
# Alternating-series control constants
C_m <- function(m, t, theta){
  
  flag <- TRUE
  
  i <- 0
  
  while (flag) {
    if (i+m == 0) {
      if( b(1, 0, t, theta) < b(0, 0, t, theta) ){ #old, but numerical problems
        flag <- FALSE
        return(i)
      }
      i <- i + 1
    }
    
    if (i+m > 0) {
      if( b_frac(i + m, m, t, theta) < 1 ){
        flag <- FALSE
        return(i)
      }
      i <- i + 1
    }
    
  }
  
}
  
# Alternating-series evaluation
S <- function(K, M = length(K)-1, t, theta, B_km){
  
  S_plus <- S_minus <- 0
  
  my_max <- max( 2*K + (1:length(K)) + 2 )
  
  if( ncol(B_km) < length(K) ){
    B_km <- cbind( B_km, matrix(data = NA, nrow = nrow(B_km), ncol = length(K) - ncol(B_km) ) )
  }
  
  if( nrow(B_km) < my_max){
    B_km <- rbind(B_km, matrix(data = NA, nrow = my_max - nrow(B_km), ncol = ncol(B_km) ) )
  }
  
  for(m in 0:M){
    
    #print(m)
    
    temp <- 0
    
    for (i in 0:( 2*K[m+1] ) ) {
      
      if( !is.na( B_km[m+i+1, m+1] ) ){
        temp <- temp + B_km[m+i+1, m+1]
      }
      else{
        temp <- temp + (-1)^i * b(m + i, m, t, theta)
        #print(c( 2*K[m+1] , m+i+1, m+1))
        B_km[m+i+1, m+1] <- (-1)^i * b(m + i, m, t, theta)
      }
      
    }
    #print(paste('temp',temp))
    
    S_plus <- S_plus + temp
    #print(paste('S_plus',S_plus))
    S_minus <- S_minus + temp
    
    #print(c(m + 2*K[m+1] + 2, m+1))
    #print(paste('DIM', dim(B_km)))
    
    if( !is.na( B_km[m + 2*K[m+1] + 2, m+1 ] ) ){
      S_minus <- S_minus + B_km[m + 2*K[m+1] + 2, m+1 ]
    }
    else{
      S_minus <- S_minus - b(m + 2*K[m+1] + 1, m, t, theta)
      B_km[m + 2*K[m+1] + 2, m+1 ] <- -b(m + 2*K[m+1] + 1, m, t, theta)
    }
    
    #print(paste('S_minus',S_minus))
    
  }
  
  return( list('minus' = S_minus, 'plus' = S_plus, 'B_km' = B_km) )
  
}
  
# Exact ancestral sampler
Ancestral <- function(t, theta, B_km = matrix(data = NA, nrow = 2, ncol = 2), nParticles = 1 ){ #from Jenkins, Spanò 2017
  
  message('    Sampling Ancestral process...', appendLF = FALSE)
  m_vec <- c()
  
  for (i in 1:nParticles) {
    
    if ( i %% 1000 == 0) message("...", i/nParticles * 100, "%", appendLF = F)

    m <- 0
    k_0 <- 0
    K <- c(k_0)
    
    U <- runif(1)
    
    flag <- TRUE
    
    while(flag){
      
      k_m <- ceiling( C_m(m, t, theta) / 2)
      
      S_temp <- S(K, m, t, theta, B_km)
      
      B_km <- S_temp$B_km
      
      while ( S_temp$minus < U && U < S_temp$plus  ) {
        K <- K + 1
        S_temp <- S(K, m, t, theta, B_km)
        B_km <- S_temp$B_km
      }
      
      if( S_temp$minus > U){
        flag <- FALSE
        m_vec <- c(m_vec, m) #new
        #return( list('m' = m, 'B_km' = B_km) ) #old
      }
      
      if(S_temp$plus < U){
        K[m+1] <- k_m
        K <- c(K, 0)
        m <- m + 1
      }

    }
  
  }
  
  message(appendLF = TRUE)
  return( list('m' = m_vec, 'B_km' = B_km) )
  
}

# Gaussian approximation to the ancestral sampler
Ancestral_small <- function(t, theta, nParticles = 1){ #from Jenkins, Spanò 2017
  #m_vec <- c()
  
  message('    Sampling Ancestral process (Gaussian approximation)...')
  #for (i in 1:nParticles) {
    #if ( i %% 100 == 0) message("...", i/nParticles * 100, "%", appendLF = F)
    
    beta <- (theta - 1) * t / 2
    
    if(beta != 0) eta <- beta / ( exp(beta) - 1 ) else eta <- 1
    
    mu <- 2*eta / t
    
    if(beta != 0) sigma2 <- mu*(eta + beta)^2 * (1 + eta/(eta+beta) - 2*eta) / beta^2 else sigma2 <- 2/(3*t)
    
    m <- rnorm(nParticles, mean = mu, sd = sqrt(sigma2) )
    
    m <- sapply(m, function(el) if (el < 0) 0 else round(el) )
    #if(m < 0) m <- 0
    
    #m_vec <- c(m_vec, round(m))
  #}

  return(m)
  
}
  
# Petrov diffusion simulation
Petrov <- function(alpha, theta, steps, delta, epsilon, N = c(10, 20, 30), show_top = Inf, start = NA, verbose = FALSE ) {
  #to make things easier let's assume times are equally spaced by delta
  #N is the fixed size of observed partition
  #m is the size, determined by the Ancestral, of omega in the formula of the transition function recalled above.
  #the smaller the delta -> the larger m -> the closer the states at t and t+delta aka continuity of the diffusion
  
  x <- pd(alpha, theta, epsilon, show_top = Inf)$weights
  
  gamma <- lapply( N, function(k) paintbox(k, x) )
  
  path <- vector(mode = 'list', length = steps + 1)
  
  path[[1]] <- list('m' = c(0), 'omega' = c(0), 'gamma' = gamma, 'x' = x)

   if (verbose) plot(x, type = 'h', main = paste('X(t0)', sep = ''), ylim = c(0,1), xlim = c(1,25), bty = 'l', xlab = 'Rank')
  
  if(delta >= 0.05) B_km = matrix(data = NA, nrow = 2, ncol = 2)
  
  for(i in 1:steps){
    
    if(delta >= 0.05) A <- Ancestral(delta, theta, B_km) else A <- Ancestral_small(delta, theta)
    
    if(delta >= 0.05) m <- A$m else m <- A
    
    if(delta >= 0.05) B_km <- A$B_km
    
    if( m == 0){
      x <- pd(alpha, theta, epsilon, show_top = Inf)$weights
      
      gamma <- lapply( N, function(k) paintbox(k, x) )
      
      path[[i+1]] <- list('m' = 0, 'omega' = c(0), 'gamma' = gamma, 'x' = x)
      
      if (verbose) plot(x, type = 'h', main = paste('X(t', i, ')', sep = ''), ylim = c(0,1), xlim = c(1,25), xlab = 'Rank')
      
      next
    }
    
    omega <- paintbox(m, x)
    
    x <- pd_posterior(alpha, theta, omega, epsilon, show_top = Inf)$weights
    
    gamma <- lapply( N, function(k) paintbox(k, x) )
    
    path[[i+1]] <- list('m' = m, 'omega' = omega, 'gamma' = gamma, 'x' = x)
    
    if (verbose) plot(x, type = 'h', main = paste('X(t', i, ')', sep = ''), ylim = c(0,1), xlim = c(1, 25), bty = 'l', xlab = 'Rank')
    
  }
  
  return(path)
  
}

Petrov_OneStep <- function(alpha, theta, delta, epsilon, X0 = NULL, show_top = Inf, min_length = 0, B_km = NULL, nParticles = 1) {
  #to make things easier let's assume times are equally spaced by delta
  #X0 must be a list of points in the simplex - each will be propagated forward in time according to Petrov's transition
  #N is the fixed size of observed partition
  #m is the size, determined by the Ancestral, of omega in the formula of the transition function recalled above.
  #the smaller the delta -> the larger m -> the closer the states at t and t+delta aka continuity of the diffusion
  
  if (is.null(X0)) X0 <- lapply(1:nParticles, function(i) pd(alpha = alpha, theta = theta, epsilon = epsilon, show_top = show_top, min_length = min_length)$weights)
  
  nParticles <- length(X0)
  
  if(delta >= 0.05 && is.null(B_km)) B_km <- matrix(data = NA, nrow = 2, ncol = 2)
  
  if(delta >= 0.05) A <- Ancestral(delta, theta, B_km, nParticles = nParticles) else A <- Ancestral_small(delta, theta, nParticles = nParticles)
  
  if(delta >= 0.05) m_vec <- A$m else m_vec <- A
  
  if(delta >= 0.05) B_km <- A$B_km
  
  message('    Propagating...')
  empirical <- pblapply(1:nParticles, function(i) {
    omega <- paintbox(m_vec[i], X0[[i]])
    x <- pd_posterior(alpha, theta, omega, epsilon, show_top = show_top, min_length = min_length)$weights
    return(x)
  })
  
  #if( m == 0){
    #x <- pd(alpha, theta, epsilon, show_top = show_top, min_length = min_length)$weights
    #return(list("X" = x) )
  #}

  if (delta >= 0.05) return(list("empirical" = empirical, "B_km" = B_km) ) else return( list("empirical" = empirical) )
  
}

#--------------------------------------------------
# Combinatorics on partitions

isLower <- function(omega, eta){
  
  l_omega <- length(omega)
  l_eta <- length(eta)
  s_omega <- sum(omega)
  s_eta <- sum(eta)
  
  if( l_omega > l_eta || s_omega > s_eta ){
    return(FALSE)
  }
  
  omega <- c(omega, rep(0, l_eta - l_omega))
  
  diff <- eta - omega
  
  if(prod(diff >= 0) == 1){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
  
}
  
# Allelic-count representation of partitions
a <- function( partition , k = sum(partition) ){
  if( length(partition) == 0) return(c())
  
  t <- table(partition)
  allelic <- rep(0, k)
  
  for(i in unique(partition)){
    allelic[i] <- t[as.character(i)]
  }
  return(allelic)
}
  
AllelicToPartition <- function( allelic ){
  partition <- c()
  
  if(length(allelic) == 0) return(c())
  
  for(i in length(allelic):1){
    if(allelic[i] > 0){
      partition <- c(partition, rep(i, allelic[i]))
    }
  }
  
  return(partition)
}

UniqueSub <- function( partition, n, as.allelic = F){
  #print(partition)
  #print(n)
  #if(length(partition) == 1 && n == 1) return(partition[1])
  
  tab <- table(partition)
  unique_sub <- comboGeneral(names(tab), m = n, freqs = tab)
  unique_sub <- lapply(1:nrow(unique_sub), function(i) sort(as.numeric(unique_sub[i, ]), decreasing = T) )
  #print(unique_sub)
  #unique_sub <- apply( unique_sub, 1, function(el) sort(as.numeric( el ), decreasing = T) )

  if( as.allelic ) unique_sub <- lapply(unique_sub, function(el) a(el, k = sum(partition)) )
  
  return( unique_sub )
}

C <- function(partition) {
  a_partition <- a(partition)
  num <- c(sum(partition):1, rep(1, length(a_partition)) )
  den_1 <- unlist(sapply(partition, function(el) el:1 ))
  den_2 <- factorial(a_partition)
  den <- c(den_1, den_2)
  return( prod(num / den) )
}

MultiComb <- function( allelic, h ){
  return( prod( CombN( allelic, h) ) )
}
  
#--------------------------------------------------
# Coagulations

CoagulateBase <- function( omega, gamma, omega2 = c(), gamma2 = c(), extra_coeff = 1 ){
  
  library(collections)
  
  l_o <- length(omega)
  l_g <- length(gamma)
  a_g <- a(gamma)
  
  coeff <- prod( factorial( a_g ) ) * extra_coeff
  
  if( l_o != l_g ) return(FALSE)
  
  coag_set <- dict()
  coag_set$set(sort( c(omega + gamma, omega2, gamma2), decreasing = T), list('kappa' = coeff ) )
  
  gamma_tab <- table(gamma)
  gamma_permutations <- permuteGeneral( names(gamma_tab), freqs = gamma_tab )
  gamma_permutations <- lapply(1:nrow(gamma_permutations), function(i) as.numeric(gamma_permutations[i, ]) )
  gamma_permutations <- do.call(rbind, gamma_permutations)
  
  if( nrow(gamma_permutations) > 1){
    for(i in 1:(nrow(gamma_permutations)-1) ){
      
      mu <- sort( c(omega + gamma_permutations[i,], omega2, gamma2), decreasing = T )
      
      kappa_old <- coag_set$get(mu, default = list('kappa' = 0))$kappa
      kappa_new <- kappa_old + coeff
      
      coag_set$set(mu, list('kappa' =  kappa_new ) )
      
    }
  }
  
  return( coag_set )
}
 
Coagulate_d <- function(omega, gamma, d, alpha, theta){
  # takes d parts from each partition and sum them
  
  coag_set <- dict()
  
  if( prod(omega == 0) == 1 ){
    coag_set$set(gamma, list('kappa' = 1))
    return(coag_set)
  }
  
  if( prod(gamma == 0) == 1){
    coag_set$set(omega, list('kappa' = 1))
    return(coag_set)
  }
  
  l_o <- length(omega)
  l_g <- length(gamma)
  
  if(l_g > l_o){
    temp <- omega
    omega <- gamma
    l_o <- l_g
    gamma <- temp
    l_g <- length(temp)
  }
  
  a_omega <- a(omega)
  a_gamma <- a(gamma)
  
  if( d > l_g) return(FALSE)
  
  
  if( d == l_o && l_o == l_g){
    coag_set <- CoagulateBase( omega, gamma )
    
    for( mu in coag_set$keys() ){
      k <- coag_set$get(mu)$kappa
      coag_set$set(mu, list('kappa' = k, 'prob' = k*PSF(alpha, theta, mu, normalized = F)))
    }
    
    return( coag_set )
  }
  
  
  if(d == 0){
    mu <- sort( c(omega, gamma), decreasing = T )
    coag_set$set(mu, list('kappa' = 1, 'prob' = PSF(alpha, theta, mu, normalized = F)) )
    return( coag_set )
  }
  
  
  if( d == l_g ){
    
    a_omega_leftout_set <- UniqueSub(omega, l_o - d, as.allelic = T)
    
    for( a_omega_leftout in a_omega_leftout_set ){
      omega_new <- AllelicToPartition( a_omega - a_omega_leftout )
      omega_out <- AllelicToPartition(a_omega_leftout)
      
      mu_base <- CoagulateBase(omega_new, gamma, omega_out, gamma2 = c(), MultiComb(a_omega, a_omega_leftout)  )
      
      for(mu in mu_base$keys() ){
        kappa_old <- coag_set$get(as.numeric(mu), default = list('kappa' = 0))$kappa
        kappa_new <- kappa_old + mu_base$get(mu)$kappa
        
        prob_old <- coag_set$get(as.numeric(mu), default = list('prob' = 0))$prob
        #prob_new <- prob_old + kappa_new*PSF(alpha, theta, mu, normalized = F) #mistake
        prob_new <- prob_old + mu_base$get(mu)$kappa*PSF(alpha, theta, mu, normalized = F) #correct
        
        coag_set$set(as.numeric(mu), list('kappa' = kappa_new, 'prob' = prob_new ))
        
      }
      
    }
    
    
    #PSF_coag <- lapply(coag_set$as_list(), function(x) x$kappa * PSF(alpha, theta, x, normalized = F))
    #sum_PSF_coag <- sum(unlist(PSF_coag))
    #PSF_coag <- lapply(PSF_coag, function(x) x / sum_PSF_coag )
    
    return( coag_set )
    
  }
  
  
  if( d < l_g){
    a_omega_leftout_set <- UniqueSub(omega, l_o - d, as.allelic = T)
    a_gamma_leftout_set <- UniqueSub(gamma, l_g - d, as.allelic = T)
    
    for( a_omega_leftout in a_omega_leftout_set ){
      
      omega_new <- AllelicToPartition( a_omega - a_omega_leftout )
      omega_out <- AllelicToPartition(a_omega_leftout)
      
      for( a_gamma_leftout in a_gamma_leftout_set ){
        gamma_new <- AllelicToPartition( a_gamma - a_gamma_leftout )
        gamma_out <- AllelicToPartition( a_gamma_leftout )
        
        mu_base <- CoagulateBase(omega_new, gamma_new, omega_out, gamma_out, MultiComb(a_omega, a_omega_leftout)*MultiComb(a_gamma, a_gamma_leftout)  )
        
        for(mu in mu_base$keys() ){
          kappa_old <- coag_set$get(as.numeric(mu), default = list('kappa' = 0))$kappa
          kappa_new <- kappa_old + mu_base$get(mu)$kappa
          
          prob_old <- coag_set$get(as.numeric(mu), default = list('prob' = 0))$prob
          #prob_new <- prob_old + kappa_new*PSF(alpha, theta, mu, normalized = F) #mistake
          prob_new <- prob_old + mu_base$get(mu)$kappa*PSF(alpha, theta, mu, normalized = F) #correct
          
          coag_set$set(as.numeric(mu), list('kappa' = kappa_new, 'prob' = prob_new))
          
        }
        
        
      }
      
    }
    
    #PSF_coag <- lapply(coag_set$keys(), function(x) coag_set$get(x)$kappa * PSF(alpha, theta, x, normalized = F))
    #sum_PSF_coag <- sum(unlist(PSF_coag))
    #PSF_coag <- lapply(PSF_coag, function(x) x / sum_PSF_coag )
    
    return( coag_set )
    
  }
  
}
 
CoagulationSet <- function(omega, gamma, alpha, theta, verbose = T){
  # Coagulate_d for d = 0, 1, 2, ...
  
  coag_set <- dict()
  
  l_o <- length(omega)
  l_g <- length(gamma)
  l <- min(l_o, l_g)
  
  if(verbose) message("    Building the set of coagulations...")
  
  S <- 0
  
  if( prod(omega == 0) == 1 ){
    coag_set$set(gamma, list('kappa' = 1, 'prob' = 1))
    return(coag_set)
  }
  
  if( prod(gamma == 0) == 1){
    coag_set$set(omega, list('kappa' = 1, 'prob' = 1))
    return(coag_set)
  }

  temp <- mclapply(0:l, function(d) Coagulate_d(omega, gamma, d = d, alpha, theta), mc.cores = detectCores() - 1 ) #use pbmclapply
  
  for( d in 1:(l+1) ){
    current <- temp[[d]]
    S_d <- unlist( lapply( current$values(), function(x) x$prob ) )
    S <- S + sum(S_d)
    coag_set$update( current )
  }
  
  temp <- NULL
  
  urn_prob <- c()
  kappas <- c()
  my_labels <- c()
  
  for(mu in coag_set$keys()){
    k <- coag_set$get(mu)$kappa
    p <- coag_set$get(mu)$prob / S
    
    urn_prob <- c(urn_prob, p)
    kappas <- c(kappas, k)
    my_labels <- c(my_labels, toString(mu))
    
    coag_set$set(mu, list('kappa' = k, 'prob' = p ))
  }
  
  if(verbose){
    plot( kappas, type = 'h', xaxt = 'n', main = paste('Kappa over Coag(', toString(omega), ';', toString(gamma), ')' ), xlab = '', ylab = 'Kappa')
    axis(1, at = 1:coag_set$size(), labels = my_labels, las = 2, cex.axis = .5)
    
    PSF_distribution <- unlist( lapply( coag_set$keys(), function(x) PSF(alpha, theta, x, normalized = FALSE)) )
    
    plot( PSF_distribution, type = 'h', ylim = c(0,max(PSF_distribution)), xaxt = 'n',
          main = paste('Unnormalised PSF over Coag(', toString(omega), ';', toString(gamma),')' ),
          sub = paste('Alpha =', alpha, 'theta=', theta) , xlab = '', ylab = 'Prob.')
    axis(1, at = 1:coag_set$size(), labels = my_labels, las = 2, cex.axis = .5)
    
    plot( urn_prob, type = 'h', ylim = c(0,max(urn_prob)), xaxt = 'n',
          main = paste('Unlabeled Pòlya Urn distribution over Coag(', toString(omega), ';', toString(gamma),')' ),
          sub = paste('Alpha =', alpha, 'theta=', theta),
          xlab = '', ylab = 'Prob.')
    axis(1, at = 1:coag_set$size(), labels = my_labels, las = 2, cex.axis = .5)
  }
  
  return(coag_set) #the length of the returned list is the cardinality of the coagulation set
}

CoagulationSetUnnormalized <- function(omega, gamma, alpha, theta, verbose = T){
  # Coagulate_d for d = 0, 1, 2, ...
  
  coag_set <- dict()
  
  l_o <- length(omega)
  l_g <- length(gamma)
  l <- min(l_o, l_g)
  
  #message("    Building the set of coagulations..." )
  
  if( prod(omega == 0) == 1 ){
    coag_set$set(gamma, list('kappa' = 1, 'prob' = PSF(alpha, theta, gamma, normalized = FALSE)))
    return(coag_set)
  }
  
  if( prod(gamma == 0) == 1){
    coag_set$set(omega, list('kappa' = 1, 'prob' = PSF(alpha, theta, omega, normalized = FALSE)))
    return(coag_set)
  }
  
  temp <- mclapply(0:l, function(d) Coagulate_d(omega, gamma, d, alpha, theta), mc.cores = detectCores() - 1 ) #use pbmclapply to print the progress bar
  
  for( d in 1:(l+1) ) coag_set$update( temp[[d]] )
  
  temp <- NULL
  
  if(verbose){
    urn_prob <- sapply(coag_set$values(), function(el) el$prob )
    kappas <- sapply(coag_set$values(), function(el) el$kappa )
    my_labels <- sapply(coag_set$keys(), function(el) toString(el) )
    
    plot( kappas, type = 'h', xaxt = 'n', main = paste('Kappa over Coag(', toString(omega), ';', toString(gamma), ')' ), xlab = '', ylab = 'Kappa')
    axis(1, at = 1:coag_set$size(), labels = my_labels, las = 2, cex.axis = .5)
    
    PSF_distribution <- unlist( lapply( coag_set$keys(), function(x) PSF(alpha, theta, x, normalized = FALSE)) )
    
    plot( PSF_distribution, type = 'h', ylim = c(0,max(PSF_distribution)), xaxt = 'n',
          main = paste('Unnormalised PSF over Coag(', toString(omega), ';', toString(gamma),')' ),
          sub = paste('Alpha =', alpha, 'theta=', theta) , xlab = '', ylab = 'Prob.')
    axis(1, at = 1:coag_set$size(), labels = my_labels, las = 2, cex.axis = .5)
    
    plot( urn_prob, type = 'h', ylim = c(0,max(urn_prob)), xaxt = 'n',
          main = paste('Unlabeled Pòlya Urn distribution over Coag(', toString(omega), ';', toString(gamma),')' ),
          sub = paste('Alpha =', alpha, 'theta=', theta),
          xlab = '', ylab = 'Prob.')
    axis(1, at = 1:coag_set$size(), labels = my_labels, las = 2, cex.axis = .5)
  }
  
  return(coag_set) #the length of the returned list is the cardinality of the coagulation set
}

#--------------------------------------------------
# Update step in the filter

Update <- function( active_nodes = dict( keys = 0, items = 1 ) , gamma, alpha, theta, pruning = 0, verbose = FALSE){
  
  coag_set <- dict()
  
  for(omega in active_nodes$keys() ){
    
    temp <- CoagulationSet(omega, gamma, alpha, theta, verbose = F)
    prob_omega <- active_nodes$get(omega)$prob
    kappa_omega <- active_nodes$get(omega)$kappa
    
    for(mu in temp$keys() ){
      
      prob_old <- coag_set$get(as.numeric(mu), default = list("kappa" = 0, "prob" = 0) )$prob
      prob_new <- prob_old + temp$get(mu)$prob * prob_omega
      
      kappa_old <- coag_set$get(as.numeric(mu), default = list("kappa" = 0, "prob" = 0) )$kappa
      kappa_new <- kappa_old + temp$get(mu)$kappa * kappa_omega
      
      coag_set$set(as.numeric(mu), list("kappa" = kappa_new, "prob" = prob_new) )
      
    }
    
  }
  
  if(pruning > 0 && pruning < 1){ #pruning is a prob
    S <- 1
    for( mu in coag_set$keys() ){
      if( coag_set$get(mu) < pruning ){
        S <- S - coag_set$get(mu)
        coag_set$pop(mu)
      }
    }
    
    for(mu in coag_set$keys() ){
      prob_old <- coag_set$get(mu)
      prob_new <- prob_old / S
      
      coag_set$set(mu, prob_new)
    }
  } 
  
  if(pruning >= 1 && pruning < coag_set$size() ){ #pruning is an integer M and take the first largest M nodes
    S <- 0
    top <- order(unlist( coag_set$values() ), decreasing = T)[1:pruning]
    coag_set_pruned <- dict()
    for(i in top){
      current_node <- coag_set$keys()[[i]]
      current_prob <- coag_set$get(coag_set$keys()[[i]])
      S <- S + current_prob 
      coag_set_pruned$set(current_node, current_prob )
    }
    
    coag_set <- coag_set_pruned
    
    for(mu in coag_set$keys() ){
      prob_old <- coag_set$get(mu)
      prob_new <- prob_old / S
      
      coag_set$set(mu, prob_new)
    }
    
  } 
  
  if(verbose){
    
    card_pre_pruning <- coag_set$size()
    print(paste('     Nodi attivi update pre pruning:', card_pre_pruning))
    
    my_labels <- unlist( lapply(coag_set$keys(), function(x) toString(x) ) )
    probs <- unlist( coag_set$values() )
    
    big_probs_id <- which(probs > 0) #which(probs >= 0.025)
    probs <- probs[big_probs_id]
    my_labels <- my_labels[big_probs_id]
    
    plot( probs, type = 'h', xaxt = 'n',
          main = paste('Update - gamma=', toString(gamma) ),
          sub = paste('alpha=', alpha, 'theta=', theta),
          xlab = '', ylab = 'Prob.')
    axis(1, at = 1:length(probs), labels = my_labels, las = 2, cex.axis = .5)
  }
  
  return(coag_set)
  
}

UpdateUnnormalized <- function( active_nodes = dict( keys = 0, items = 1 ) , gamma, alpha, theta, pruning = 0, verbose = FALSE){
  
  coag_set <- dict()
  
  for(omega in active_nodes$keys() ){
    
    temp <- CoagulationSetUnnormalized(omega, gamma, alpha, theta, verbose = F)
    prob_omega <- active_nodes$get(omega)$prob
    kappa_omega <- active_nodes$get(omega)$kappa
    
    for(mu in temp$keys() ){
      
      prob_old <- coag_set$get(as.numeric(mu), default = list("kappa" = 0, "prob" = 0) )$prob
      prob_new <- prob_old + temp$get(mu)$prob * prob_omega
      
      kappa_old <- coag_set$get(as.numeric(mu), default = list("kappa" = 0, "prob" = 0) )$kappa
      kappa_new <- kappa_old + temp$get(mu)$kappa * kappa_omega
      
      coag_set$set(as.numeric(mu), list("kappa" = kappa_new, "prob" = prob_new) )
      
    }
    
  }
  
  if(pruning > 0 && pruning < 1){ #pruning is a prob
    S <- 1
    for( mu in coag_set$keys() ){
      if( coag_set$get(mu) < pruning ){
        S <- S - coag_set$get(mu)
        coag_set$pop(mu)
      }
    }
    
    for(mu in coag_set$keys() ){
      prob_old <- coag_set$get(mu)
      prob_new <- prob_old / S
      
      coag_set$set(mu, prob_new)
    }
  } 
  
  if(pruning >= 1 && pruning < coag_set$size() ){ #pruning is an integer M and take the first largest M nodes
    S <- 0
    top <- order(unlist( coag_set$values() ), decreasing = T)[1:pruning]
    coag_set_pruned <- dict()
    for(i in top){
      current_node <- coag_set$keys()[[i]]
      current_prob <- coag_set$get(coag_set$keys()[[i]])
      S <- S + current_prob 
      coag_set_pruned$set(current_node, current_prob )
    }
    
    coag_set <- coag_set_pruned
    
    for(mu in coag_set$keys() ){
      prob_old <- coag_set$get(mu)
      prob_new <- prob_old / S
      
      coag_set$set(mu, prob_new)
    }
    
  } 
  
  if(verbose){
    
    card_pre_pruning <- coag_set$size()
    print(paste('     Nodi attivi update pre pruning:', card_pre_pruning))
    
    my_labels <- unlist( lapply(coag_set$keys(), function(x) toString(x) ) )
    probs <- unlist( coag_set$values() )
    
    big_probs_id <- which(probs > 0) #which(probs >= 0.025)
    probs <- probs[big_probs_id]
    my_labels <- my_labels[big_probs_id]
    
    plot( probs, type = 'h', xaxt = 'n',
          main = paste('Update - gamma=', toString(gamma) ),
          sub = paste('alpha=', alpha, 'theta=', theta),
          xlab = '', ylab = 'Prob.')
    axis(1, at = 1:length(probs), labels = my_labels, las = 2, cex.axis = .5)
  }
  
  return(coag_set)
  
}

UpdateUnnormalized2 <- function( active_nodes = dict( keys = 0, items = 1 ) , gamma, alpha, theta, epsilon = 0.005, N = 10^4, show_top = 3, pruning = 0, verbose = FALSE){
  
  coag_set <- dict()
  
  #message('    Update...')
  lapply(active_nodes$keys(), function(omega) { #use pblapply
    temp <- CoagulationSetUnnormalized(omega, gamma, alpha, theta, verbose = F)
    prob_omega <- active_nodes$get(omega)
    
      for(mu in temp$keys() ){
        
        prob_old <- coag_set$get(as.numeric(mu), default = 0 )
        prob_new <- prob_old + temp$get(mu)$prob * prob_omega / PSF(alpha, theta, omega, normalized = FALSE)
        
        coag_set$set(as.numeric(mu), prob_new )
        
      }
    
    })
  
  #message(paste('    Active nodes after update (pre-pruning):', coag_set$size()) )
  
  S <- sum(unlist(coag_set$values()))
  
  coag_set <- dict(keys = coag_set$keys(), items = as.list( unlist(coag_set$values()) / S ))
  
  #eliminate this if you want to sample from the pruned mixture (fix accordingly the function Particle_filter)
  if (N > 0) {
    #message(paste('    Drawing from the filtering distribution (before pruning)...'))
    
    #message('    Unlist...', appendLF = FALSE)
    coag_set_distribution <- unlist( coag_set$values() )
    #message('done!', appendLF = T)
    #message('    Sampling...', appendLF = FALSE)
    coag_set_index <- sample(1:coag_set$size(), N, replace = TRUE, prob = coag_set_distribution)
    #message('done!', appendLF = T)
    #print(sort(as.numeric(table(coag_set_index)), decreasing = T))
    #message('    Lapply...', appendLF = FALSE)
    #mus <- pblapply(coag_set_index, function(i) coag_set$keys()[[i]] )
    mus <- coag_set$keys()[coag_set_index]
    #message('done!', appendLF = T)
    #message('    Posterior...', appendLF = FALSE)
    filter <- t( sapply(mus, function(mu) pd_posterior(alpha, theta, mu, epsilon, show_top = show_top )$weights ) ) #use pbsapply
    #message('done!', appendLF = T)
    
    #message('    % of unique sampled nodes: ', round(length(unique(mus)) / coag_set$size(), 2) )
  }
  
  if (FALSE) {
    coag_lengths <- sapply(coag_set$keys(), function(el) length(el))
    coag_lengths <- data.frame("lengths" = coag_lengths, "probs" = unlist(coag_set$values()) )
    coag_lengths <- coag_lengths %>% group_by(lengths) %>% summarise_at(vars(probs), list(probs = sum))
    my_order <- order(coag_lengths$probs, decreasing = TRUE)
    print(coag_lengths[my_order, ])
  }
  
  #pruning is a prob: keep only the nodes having largest probability with cumulative probability at least [pruning]
  if(pruning > 0 && pruning < 1){
    #message('    Pruning...', appendLF = FALSE)
    top_index <- order(unlist( coag_set$values() ), decreasing = T)
    top_nodes_prob <- unlist(coag_set$values())[top_index]
    cumulative <- cumsum(top_nodes_prob)
    cutoff <- sum(cumulative < pruning) + 1
    top_nodes <- lapply(top_index[1:cutoff], function(i) coag_set$keys()[[i]] )
    top_nodes_prob <- top_nodes_prob[1:cutoff]
    coag_set <- dict(keys = top_nodes, items = as.list(top_nodes_prob / sum(top_nodes_prob) ))
    
    message(paste(' active nodes post-pruning:', coag_set$size()))
  } 
  
  #pruning is an integer M and take the first largest M nodes
  if(pruning >= 1 && pruning < coag_set$size() ){
    #message('    Pruning...')
    
    top_index <- order(unlist( coag_set$values() ), decreasing = T)[1:pruning]
    top_nodes <- lapply(top_index, function(i) coag_set$keys()[[i]] )
    top_nodes_prob <- lapply(top_index, function(i) coag_set$values()[[i]] )
    S <- sum(unlist(top_nodes_prob))
    top_nodes_prob <- lapply(top_nodes_prob, function(el) el / S)
    coag_set <- dict(items = top_nodes_prob, keys = top_nodes)
  }
  
  if(verbose){
    
    my_labels <- unlist( lapply(coag_set$keys(), function(x) toString(x) ) )
    probs <- unlist( coag_set$values() )
    my_order <- order(probs, decreasing = TRUE)
    probs <- probs[my_order]
    my_labels <- my_labels[my_order]
    
    par(mfrow = c(1, 1), mar = c(6, 4, 4, 1))
    plot( probs, type = 'h', xaxt = 'n',
          main = paste('Update - gamma=', toString(gamma) ),
          xlab = '', ylab = 'Prob.', bty = 'l' )
    axis(1, at = 1:length(probs), labels = my_labels, las = 2, cex.axis = .33)
  }
  
  if (N > 0) return(list(coag_set, filter)) else return(list(coag_set, 'lik' = S))
  
}

multiUpdate <- function( partition.list, alpha, theta, pruning = 0, verbose = FALSE, likelihood = FALSE){
  
  coag_set <- ordered_dict()
  coag_set$set(partition.list[[1]], list("kappa" = 1, "prob" = 1))
  L <- length(partition.list)
  
  for (j in 2:L) {
    gamma_j <- partition.list[[j]]
    coag_set <- Update(coag_set, gamma_j, alpha, theta, verbose = FALSE)
    #if (j < L) coag_set <- Update(coag_set, gamma_j, alpha, theta, verbose = FALSE)
    #if (j == L) coag_set <- UpdateUnnormalized(coag_set, gamma_j, alpha, theta, verbose = FALSE)
    #if ( (100 * j / L) %% 1 == 0 ) message(paste0(j / L * 100, "%..."), appendLF = FALSE)
    message(paste0(j / L * 100, "%..."), appendLF = FALSE)
  }
  
  for(mu in coag_set$keys()){
    k <- coag_set$get(mu)$kappa
    p <- k * PSF(alpha, theta, mu, normalized = FALSE)
    coag_set$set(mu, list("kappa" = k, "prob" = p))
  }
  
  S <- sum(sapply(coag_set$values(), function(val) val$prob ))
  
  for(mu in coag_set$keys()){
    k <- coag_set$get(mu)$kappa
    p <- coag_set$get(mu)$prob / S
    coag_set$set(mu, list("kappa" = k, "prob" = p))
  }
  
  if(pruning > 0 && pruning < 1){ #pruning is a prob
    S <- 1
    for( mu in coag_set$keys() ){
      if( coag_set$get(mu) < pruning ){
        S <- S - coag_set$get(mu)
        coag_set$pop(mu)
      }
    }
    
    for(mu in coag_set$keys() ){
      prob_old <- coag_set$get(mu)
      prob_new <- prob_old / S
      
      coag_set$set(mu, prob_new)
    }
  } 
  
  if(pruning >= 1 && pruning < coag_set$size() ){ #pruning is an integer M and take the first largest M nodes
    S <- 0
    top <- order(unlist( coag_set$values() ), decreasing = T)[1:pruning]
    coag_set_pruned <- dict()
    for(i in top){
      current_node <- coag_set$keys()[[i]]
      current_prob <- coag_set$get(coag_set$keys()[[i]])
      S <- S + current_prob 
      coag_set_pruned$set(current_node, current_prob )
    }
    
    coag_set <- coag_set_pruned
    
    for(mu in coag_set$keys() ){
      prob_old <- coag_set$get(mu)
      prob_new <- prob_old / S
      
      coag_set$set(mu, prob_new)
    }
    
  } 
  
  if(verbose){
    #extract the coagulations (support of the target posterior)
    label_lengths <- sapply(coag_set$keys(), function(el) length(el))
    my_labels <- sapply( coag_set$keys(), function(el) toString(el))
    L <- length(my_labels)
    
    #sort the support from the longest to the shortes partition
    my_labels <- my_labels[order(label_lengths, decreasing = TRUE)]
    
    #extract the target kappa*eppf (it's already normalized...)
    prob <- sapply( coag_set$values(), function(val) val$prob)
    prob <- prob[order(label_lengths, decreasing = TRUE)]
    names(prob) <- my_labels
    
    #extract the kappas
    kappa <- sapply( coag_set$values(), function(val) val$kappa)
    kappa_norm <- kappa / sum(kappa)
    kappa <- kappa[order(label_lengths, decreasing = TRUE)]
    names(kappa) <- my_labels
    
    #PSF a priori
    prior <- sapply(coag_set$keys(), function(mu) PSF(alpha = alpha, theta = theta, partition = mu, normalized = TRUE) )
    prior <- prior[order(label_lengths, decreasing = TRUE)]
    names(prior) <- my_labels
    
    par(mfrow = c(3, 1), mar = c(.5, 4, 3, 1))
    
    plot( kappa, type = 'h', xaxt = 'n',
          main = 'multiUpdate',
          xlab = '', ylab = 'Kappa')
    
    par(mar = c(.5, 4, 1, 1))
    
    plot( prior, type = 'h', xaxt = 'n',
          main = '',
          xlab = '', ylab = 'Prior')
    
    par(mar = c(6, 4, 1, 1))
    plot( prob, type = 'h', xaxt = 'n',
          main = '',
          sub = paste('alpha=', alpha, 'theta=', theta),
          xlab = '', ylab = 'Posterior')
    axis(1, at = 1:length(prob), labels = my_labels, las = 2, cex.axis = .5)
  }
  
  #return(list(coag_set, S))
  if(likelihood) return(S) else return( coag_set )
  
}

#--------------------------------------------------

Polya_urn <- function(alpha, theta, omega, gamma){
  
  #if( current > 0 && ( current / nParticles * 100 ) %% 5 == 0 ) print( paste('    Polya Urn status:', current / nParticles * 100, '%' ) )
  
  l_omega <- length(omega)
  
  s_gamma <- sum(gamma)
  s_omega <- sum(omega)
  
  if(s_omega > 0) new_color <- l_omega + 1 else new_color <- 1
  
  if( s_omega > 0) urn <- omega else urn <- c()
  
  control_urn <- c()
  
  if( s_omega > 0) probs <- urn - alpha else probs <- c()
  
  urn <- c(urn, 1)
  probs <- c(probs, theta + (new_color - 1)*alpha)
  
  count <- 0
  #iterations <- 0
  flag <- TRUE
  
  while( flag ){
    
    #iterations <- iterations + 1
    
    one_sample <- sample(1:new_color, 1, replace = TRUE, prob = probs)
    control_urn <- c(control_urn, one_sample)
    
    if( !isLower( sort( as.numeric( table(control_urn) ), decreasing = TRUE ) , gamma ) ){
      
      flag <- TRUE
      
      if(s_omega > 0) new_color <- l_omega + 1 else new_color <- 1
      
      if( s_omega > 0) urn <- omega else urn <- c()
      
      control_urn <- c()
      
      if( s_omega > 0) probs <- urn - alpha else probs <- c()
      
      urn <- c(urn, 1)
      probs <- c(probs, theta + (new_color - 1)*alpha)
      
      count <- 0
      
      next
    }
    
    count <- count + 1
    
    if( one_sample == new_color ){
      probs[new_color] <- 1 - alpha
      new_color <- new_color + 1
      urn <- c(urn, 1)
      probs <- c(probs, theta + (new_color - 1)*alpha)
    }
    else{
      urn[one_sample] <- urn[one_sample] + 1
      probs[one_sample] <- probs[one_sample] + 1
    }
    
    if( count == s_gamma ){
      control_urn <- sort( as.numeric( table(control_urn) ), decreasing = TRUE )
      if( length(control_urn) == length(gamma) ) check <- prod( control_urn == gamma ) else check <- 0
    }
    
    if( count == s_gamma && check == 1  ){
      flag <- FALSE
    }
    
    if( count == s_gamma && check == 0 ){
      
      flag <- TRUE
      
      if(s_omega > 0) new_color <- l_omega + 1 else new_color <- 1
      
      if( s_omega > 0) urn <- omega else urn <- c()
      
      control_urn <- c()
      
      if( s_omega > 0) probs <- unlist( lapply(urn, function(x) if(x > 0) x - alpha ) ) else probs <- c()
      
      urn <- c(urn, 1)
      probs <- c(probs, theta + (new_color - 1)*alpha)
      
      count <- 0
    }
    
  }
  
  #if( floor(current / particles * 100) %% 20 == 0 && floor(current / particles * 100) != floor(old_current / particles * 100) ) message( paste0( round( current/particles * 100, 2), '%' ), appendLF = FALSE )
  
  urn <- urn[ -new_color ]
  urn <- sort(urn, decreasing = TRUE)
  #urn <- sort( as.numeric(table(urn)), decreasing = TRUE)
  
  #print(urn)
  
  return( urn )
}

Polya_urn_N <- function(alpha, theta, omega, gamma, N = 10^4) {
  message('       Polya Urn status: ', appendLF = FALSE)
  oldcurrent <- 0
  mu_list <- list()
  for(current in 1:N){
    mu_list <- c(mu_list, list(Polya_urn(alpha, theta, omega, gamma)))
    if( floor(current / N * 100) - floor(oldcurrent / N * 100) == 10  ) {
      message( paste0( round( current/N * 100), '%...' ), appendLF = FALSE )
      oldcurrent <- current
    }
    
  }
  return(mu_list)
}

Polya_urn2 <- function(alpha, theta, omega, n){
  #this samples pi, not a coagulation
  #in particular, if omega = c(), then this is the standard CRP starting from empty restaurant
  
  #if( current > 0 && ( current / nParticles * 100 ) %% 5 == 0 ) print( paste('    Polya Urn status:', current / nParticles * 100, '%' ) )
  
  l_omega <- length(omega)
  
  s_omega <- sum(omega)
  
  if(s_omega > 0) new_color <- l_omega + 1 else new_color <- 1
  
  if( s_omega > 0) urn <- omega else urn <- c()
  
  control_urn <- c()
  
  if( s_omega > 0) probs <- urn - alpha else probs <- c()
  
  urn <- c(urn, 1)
  probs <- c(probs, theta + (new_color - 1)*alpha)
  
  count <- 0
  #iterations <- 0
  flag <- TRUE
  
  while( count < n ){
    
    #iterations <- iterations + 1
    
    one_sample <- sample(1:new_color, 1, replace = TRUE, prob = probs)
    control_urn <- c(control_urn, one_sample)
    count <- count + 1
    
    if( one_sample == new_color ){
      probs[new_color] <- 1 - alpha
      new_color <- new_color + 1
      urn <- c(urn, 1)
      probs <- c(probs, theta + (new_color - 1)*alpha)
    }
    else{
      urn[one_sample] <- urn[one_sample] + 1
      probs[one_sample] <- probs[one_sample] + 1
    }
    
  }
  
  #if( floor(current / particles * 100) %% 20 == 0 && floor(current / particles * 100) != floor(old_current / particles * 100) ) message( paste0( round( current/particles * 100, 2), '%' ), appendLF = FALSE )
  
  control_urn <- sort( as.numeric( table(control_urn) ), decreasing = TRUE )
  urn <- urn[-length(urn)]
  urn <- sort(urn, decreasing = T)
  
  return( list(urn, control_urn) )
}

#--------------------------------------------------
# The dual process

a_knl <- function(k, n, l, theta){
  
  num1 <- log( theta + 2*k - 1) - sum( log( 1:l ) )
  if( k > l ) num1 <- num1 - sum( log( 1:(k - l) ) )
  
  #num2 <- 0
  
  #num2 <- -sum( log( 1:(k - l) ) )
  if(k > 1){
    for(i in 0:(k-2)){
      #num2 <- num2 + log(l + theta + i)
      num1 <- num1 + log(l + theta + i)
    }
  }
  
  #print(num1)
  #num3 <- 0
  #num3 <- sum(log( 1:n ) ) - sum( log( 1:(n-k) ) )
  num3 <- sum(log( (n-k+1):n ))
  #print(num3)
  
  if(k > 0){
    for(i in 0:(k-1)){
      num3 <- num3 - log( theta + n + i )
    }
  }
  
  #num3 <- log(factorial(n)) - log(factorial(n-k)) - num3
  
  #print(num1 + num2 + num3)
  
  return( exp(num1 + num3) )
  
  #print(num1+num2)
  
}

b_knl <- function(k, n, l, t, theta){
  
  fact1 <- a_knl(k, n, l, theta)
  
  fact2 <- exp(- ( k * (k + theta - 1) * t ) / 2 )
  
  #print(fact1*fact2)
  
  return( fact1 * fact2 )
  
}

Ancestral_transition <- function(n, l, t, theta ){
  
  if(l == 0) return(F)
  
  d <- 0
  
  for(k in l:n ){
    d <- d + (-1)^(k - l) * b_knl(k, n, l, t, theta)
  }
  
  return(d)
} #for l >= 1
 
Ancestral_transition_death <- function(n, t, theta){
  
  s <- 0
  
  for(l in 1:n){
    s <- s + Ancestral_transition(n, l, t, theta)
  }
  
  return(1-s)
  
} #i.e. when l = 0
 
rBlockCounting_single <- function(n, t, theta){
  elapsed_t <- 0
  current_n <- n
  while (elapsed_t < t) {
    Exp_t <- rexp(1, rate = current_n*(current_n + theta - 1)/2 )
    elapsed_t <- elapsed_t + Exp_t
    if (elapsed_t < t) current_n <- current_n - 1 else return(current_n)
    if (current_n == 0) return(current_n)
    if (theta == 0 && current_n == 1) return(current_n)
  }
}

rBlockCounting <- function(N, n, t, theta){
  return( sapply(1:N, function(i) rBlockCounting_single(n, t, theta) ) )
}
 
DownChain <- function(eta, t, theta){
  
  m <- rBlockCounting_single(sum(eta), t, theta)
  
  if(m == sum(eta)) return(eta)
  
  urn <- rep( 1:length(eta), eta)
  omega <- sample(urn, m)
  omega <- as.numeric( table(omega) )
  omega <- sort(omega, decreasing = T)
  
  if(length(omega) > 0) return(omega) else return(0)
  
}
 
# dual process propagation of particles (Gillespie approximation)
Particle_predict <- function( nodes = dict(items = 1, keys = c(0) ), t, theta, pruning = 0, nParticles = 10^4, verbose = FALSE, bootstrap = FALSE ){
  
  #message('    Propagating...', appendLF = FALSE)
  
  nodes_distribution <- unlist( nodes$values() )

  particles <- sample(1:nodes$size(), nParticles, replace = TRUE, prob = nodes_distribution)
  
  #etas <- lapply(particles, function(x) nodes$keys()[[x]] )
  etas <- nodes$keys()[particles]
  
  #message( paste( 'unique starting urns:', length(unique(etas)) ) )
  
  omegas <- lapply(etas, function(el) DownChain( el, t, theta) ) #use pblapply
  
  predict_set <- dict()
  for (omega in omegas) {
    old_count <- predict_set$get(omega, default = 0)
    predict_set$set(omega, old_count + 1)
  }
  
  size <- predict_set$size()
  #message(paste('    Active nodes after propagation (pre-pruning):',  size ) )
  #if(verbose) print(paste('   ', round(nParticles / predict_set$size(), 2), ' particles per node (on average)' ) )
  
  #pruning is a prob: keep the most probable nodes with cumulative prob of at least [pruning]
  if(pruning > 0 && pruning < 1){
    #message('    Pruning...', appendLF = FALSE)
    
    top_index <- order(unlist( predict_set$values() ), decreasing = T)
    top_nodes_prob <- unlist(predict_set$values())[top_index] / sum(unlist(predict_set$values()))
    cumulative <- cumsum(top_nodes_prob)
    cutoff <- sum( cumulative < pruning ) + 1
    top_nodes <- lapply(top_index[1:cutoff], function(i) predict_set$keys()[[i]] )
    top_nodes_prob <- top_nodes_prob[1:cutoff]
    predict_set <- dict(keys = top_nodes, items = as.list(top_nodes_prob / sum(top_nodes_prob) ))
    
    #message(paste(' active nodes post-pruning:', predict_set$size()))
  } 
  
  #pruning is an integer M and take the first largest M nodes
  if(pruning >= 1 && pruning < size ){ 
    #message('    Pruning...')
    
    top_index <- order(unlist( predict_set$values() ), decreasing = T)[1:pruning]
    #top_nodes <- lapply(top_index, function(i) predict_set$keys()[[i]] )
    top_nodes <- predict_set$keys()[top_index]
    #top_nodes_prob <- lapply(top_index, function(i) predict_set$values()[[i]] )
    top_nodes_prob <- predict_set$values()[top_index]
    S <- sum(unlist(top_nodes_prob))
    top_nodes_prob <- lapply(top_nodes_prob, function(el) el / S)
    predict_set <- dict(items = top_nodes_prob, keys = top_nodes)
  }
  
  if (pruning == 0) {
    check_sum <- sum( unlist( predict_set$values() ) )
    if (check_sum != nParticles) print( paste('WRONG PREDICT PRUNING!!!!!!!!!!!!', nParticles - check_sum ) )
    
    top_nodes_prob <- predict_set$values()
    S <- sum(unlist(top_nodes_prob))
    top_nodes_prob <- lapply(top_nodes_prob, function(el) el / S)
    predict_set <- dict(items = top_nodes_prob, keys = predict_set$keys() )
  }
  
  if(verbose){
    my_labels <- unlist( lapply( predict_set$keys(), function(x) toString(x) ) )
    probs <- unlist( predict_set$values() )
    my_order <- order(probs, decreasing = TRUE)
    probs <- probs[my_order]
    my_labels <- my_labels[my_order]
    
    par(mfrow = c(1, 1), mar = c(6, 4, 4, 1))
    plot( probs, type = 'h', xaxt = 'n',
          main = 'Predict',
          xlab = '', ylab = 'Prob.', bty = 'l')
    axis(1, at = 1:length(probs), labels = my_labels, las = 2, cex.axis = .33)
  }
  
  #print(predict_set$as_list())

  return( predict_set )
}

#--------------------------------------------------
# Duality-based filter

Particle_filter <- function( observation_list, alpha, theta, t, epsilon, pruning_update = 0, pruning_propagation = 0, N = 10^4, nParticles = 10^4, show_top = 5, verbose = FALSE, store_mixture_update = FALSE, store_mixture_propagation = FALSE){
  
  #N is the number of Monte Carlo replicates per observed partition
  start.time.PF <- Sys.time()
  
  nodes <- dict()
  filter <- vector(mode = 'list', length = length(observation_list))
  if(store_mixture_update == TRUE) collection_nodes_update <- vector(mode = 'list', length = length(observation_list))
  if(store_mixture_propagation == TRUE) collection_nodes_propagation <- vector(mode = 'list', length = length(observation_list) - 1)
  
  for( i in 1:length(observation_list) ){
    
    #aux.time <- Sys.time()
    
    gamma <- observation_list[[i]]
    
    message(paste("Time", i, 'out of', length(observation_list)))
    
    if(i == 1){
      nodes <- dict()
      nodes$set(gamma, nParticles)
      filter[[i]] <- t( pbsapply(1:N, function(n) pd_posterior(alpha, theta, gamma, epsilon, show_top = show_top )$weights ) )
    }
    else{
      temp <- UpdateUnnormalized2(nodes, gamma, alpha, theta, epsilon = epsilon, N = N, show_top = show_top, pruning = pruning_update, verbose = verbose)
      nodes <- temp[[1]]
      filter[[i]] <- temp[[2]]
    }
    if(store_mixture_update == TRUE) collection_nodes_update[[i]] <- nodes
    
    #partial.elapsed.time <- round(difftime(cur.time, aux.time), 3)
    #message(sprintf('Update in %.2f %s', partial.elapsed.time, units(partial.elapsed.time)))
    #elapsed.time <- round(difftime(sys.time(), start.time.PF), 3)
    #message('    Overall in ', elapsed.time, '', units(elapsed.time))
    
    if( nodes$size() == 0 && pruning > 0 && pruning < 1){
      message('Attenzione! Pruning troppo alto! Tutti i nodi sono stati eliminati!')
      return(NULL)
    }
    
    #diagnostics: the distributions of the probabilities of active nodes
    #current_probs <- sort( unlist( nodes$values() ), decreasing = T )
    #plot(current_probs, type = 'h', main = 'Update', ylim = c(0, max(current_probs)) )
    #abline(h = pruning, lty = 'dashed', col = 'red')
    
    if( i < length(observation_list) ) nodes <- Particle_predict(nodes, t, theta, pruning = pruning_propagation, nParticles = nParticles, verbose = verbose)
    if(i < length(observation_list) && store_mixture_propagation == TRUE) collection_nodes_propagation[[i]] <- nodes
    
    if( nodes$size() == 0 && pruning > 0 && pruning < 1 ){
      message('Attenzione! Pruning troppo alto! Tutti i nodi sono stati eliminati!')
      return(NULL)
    }
    #diagnostics: the distributions of the probabilities of active nodes
    #current_probs <- sort( unlist( nodes$values() ), decreasing = T )
    #plot(current_probs, type = 'h', main = 'Propagation', ylim = c(0, max(current_probs)))
    #abline(h = pruning, lty = 'dashed', col = 'red')
    
  }
  
  elapsed.time <- round(difftime(Sys.time(), start.time.PF), 3)
  message('    Duality filter runtime: ', elapsed.time, ' ', units(elapsed.time) )
  
  if(store_mixture_update == TRUE && store_mixture_propagation == TRUE) return(list(filter, collection_nodes_update, collection_nodes_propagation))
  if(store_mixture_update == TRUE && store_mixture_propagation == FALSE) return(list(filter, collection_nodes_update))
  if(store_mixture_update == FALSE && store_mixture_propagation == FALSE) return(filter)
}

Particle_filter_likelihood <- function( observation_list, alpha, theta, t, epsilon, pruning_update = 0, pruning_propagation = 0, N = 10^4, nParticles = 10^4, show_top = 5, verbose = FALSE, store_mixture = FALSE){
  
  #N is the number of Monte Carlo replicates per observed partition
  start.time.PF <- Sys.time()
  
  nodes <- dict()
  likelihood <- vector(mode = 'numeric', length = length(observation_list))
  if(store_mixture == TRUE) collection_nodes <- vector(mode = 'list', length = length(observation_list))
  
  for( i in 1:length(observation_list) ){
    
    #aux.time <- Sys.time()
    
    gamma <- observation_list[[i]]
    
    #message(paste("Time", i, 'out of', length(observation_list)))
    
    if(i == 1){
      nodes <- dict()
      nodes$set(gamma, nParticles)
      likelihood[i] <- PSF(alpha, theta, gamma, normalized = TRUE)
    }
    else{
      temp <- UpdateUnnormalized2(nodes, gamma, alpha, theta, epsilon = epsilon, N = 0, show_top = show_top, pruning = pruning_update, verbose = verbose)
      nodes <- temp[[1]]
      likelihood[i] <- temp$lik
    }
    if(store_mixture == TRUE) collection_nodes[[i]] <- nodes
    
    #partial.elapsed.time <- round(difftime(cur.time, aux.time), 3)
    #message(sprintf('Update in %.2f %s', partial.elapsed.time, units(partial.elapsed.time)))
    #elapsed.time <- round(difftime(sys.time(), start.time.PF), 3)
    #message('    Overall in ', elapsed.time, '', units(elapsed.time))
    
    if( nodes$size() == 0 && pruning > 0 && pruning < 1){
      message('Attenzione! Pruning troppo alto! Tutti i nodi sono stati eliminati!')
      return(NULL)
    }
    
    #diagnostics: the distributions of the probabilities of active nodes
    #current_probs <- sort( unlist( nodes$values() ), decreasing = T )
    #plot(current_probs, type = 'h', main = 'Update', ylim = c(0, max(current_probs)) )
    #abline(h = pruning, lty = 'dashed', col = 'red')
    
    if( i < length(observation_list) ) nodes <- Particle_predict(nodes, t, theta, pruning = pruning_propagation, nParticles = nParticles, verbose = verbose)
    
  }
  
  elapsed.time <- round(difftime(Sys.time(), start.time.PF), 3)
  message('    Duality filter runtime: ', elapsed.time, ' ', units(elapsed.time) )
  
  if(store_mixture == TRUE) return(list(likelihood, collection_nodes)) else return(likelihood)
}

Particle_filter_qentropy <- function( observation_list, alpha, theta, t, epsilon, pruning_update = 0, pruning_propagation = 0, N = 10^4, nParticles = 10^4, show_top = Inf, verbose = FALSE, store_mixture = FALSE, store_mixture_size = FALSE, q = 1){
  
  #N is the number of Monte Carlo replicates per observed partition
  start.time.PF <- Sys.time()
  
  nodes <- dict()
  filter <- vector(mode = 'list', length = length(observation_list))
  if(store_mixture) collection_nodes <- vector(mode = 'list', length = length(observation_list))
  if(store_mixture_size) collection_nodes_size <- vector(mode = 'numeric', length = length(observation_list))
  
  for( i in 1:length(observation_list) ){
    
    #aux.time <- Sys.time()
    
    gamma <- observation_list[[i]]
    
    message(paste("Time", i, 'out of', length(observation_list)))
    
    if(i == 1){
      nodes <- dict()
      nodes$set(gamma, nParticles)
      temp <- t( sapply(1:N, function(n) pd_posterior(alpha, theta, gamma, epsilon, show_top = Inf )$weights ) ) #use pbsapply
      filter[[i]] <- sapply(temp, function(x) qentropy(x / sum(x), q)  )
    }
    else{
      temp <- UpdateUnnormalized2(nodes, gamma, alpha, theta, epsilon = epsilon, N = N, show_top = Inf, pruning = pruning_update, verbose = verbose)
      nodes <- temp[[1]]
      temp <- temp[[2]]
      filter[[i]] <- sapply(temp, function(x) qentropy(x / sum(x), q = q)  )
    }
    if(store_mixture) collection_nodes[[i]] <- nodes
    if(store_mixture_size) collection_nodes_size[i] <- nodes$size()
    
    #partial.elapsed.time <- round(difftime(cur.time, aux.time), 3)
    #message(sprintf('Update in %.2f %s', partial.elapsed.time, units(partial.elapsed.time)))
    #elapsed.time <- round(difftime(sys.time(), start.time.PF), 3)
    #message('    Overall in ', elapsed.time, '', units(elapsed.time))
    
    if( nodes$size() == 0 && pruning > 0 && pruning < 1){
      message('Warning! Pruning is too high! All nodes deleted!')
      return(NULL)
    }
    
    #diagnostics: the distribution of the probabilities of active nodes
    #current_probs <- sort( unlist( nodes$values() ), decreasing = T )
    #plot(current_probs, type = 'h', main = 'Update', ylim = c(0, max(current_probs)) )
    #abline(h = pruning, lty = 'dashed', col = 'red')
    
    if( i < length(observation_list) ) nodes <- Particle_predict(nodes, t, theta, pruning = pruning_propagation, nParticles = nParticles, verbose = verbose)
    
    #diagnostics: the distribution of the probabilities of active nodes
    #current_probs <- sort( unlist( nodes$values() ), decreasing = T )
    #plot(current_probs, type = 'h', main = 'Propagation', ylim = c(0, max(current_probs)))
    #abline(h = pruning, lty = 'dashed', col = 'red')
    
  }
  
  elapsed.time <- round(difftime(Sys.time(), start.time.PF), 3)
  message('    Duality filter runtime: ', elapsed.time, ' ', units(elapsed.time) )
  
  if(store_mixture) return(list(filter, collection_nodes))
  if(store_mixture_size) return(list(filter, collection_nodes_size))
  else return(filter)
}

Particle_predictive_qentropy <- function( observation_list, alpha, theta, t, epsilon, pruning_update = 0, pruning_propagation = 0, N = 10^4, nParticles = 10^4, show_top = Inf, verbose = FALSE, store_mixture = FALSE, q = 2){
  
  #N is the number of Monte Carlo replicates per observed partition
  start.time.PF <- Sys.time()
  
  predictive <- vector(mode = 'list', length = length(observation_list) + 1)
  if(store_mixture == TRUE) collection_nodes <- vector(mode = 'list', length = length(observation_list))
  
  for( i in 1:length(observation_list) ){
    
    #aux.time <- Sys.time()
    
    message(paste("Time", i, 'out of', length(observation_list) + 1 ))
    
    if(i == 1){
      nodes <- dict( keys = 0, items = 1)
      temp <- t( pbsapply(1:N, function(n) pd(alpha, theta, epsilon, show_top = Inf )$weights ) )
      predictive[[i]] <- sapply(temp, function(x) qentropy(x / sum(x), q)  )
    }
    else{
      nodes <- Particle_predict(nodes, t, theta, pruning = 0, nParticles = nParticles, verbose = verbose)
      message(paste('    Drawing from the predictive distribution (pre-pruning)...'))
      nodes_distribution <- unlist( nodes$values() )
      nodes_index <- sample(1:nodes$size(), N, replace = TRUE, prob = nodes_distribution)
      mus <- nodes$keys()[nodes_index]
      temp <- t( pbsapply(mus, function(mu) pd_posterior(alpha, theta, mu, epsilon, show_top = show_top )$weights ) )
      message('    % of unique sampled nodes: ', round(length(unique(mus)) / nodes$size(), 2) )
      predictive[[i]] <- sapply(temp, function(x) qentropy(x / sum(x), q)  )
      
      message('    Pruning...')
      top_index <- order(unlist( nodes$values() ), decreasing = T)[1:pruning_propagation]
      top_nodes <- nodes$keys()[top_index]
      top_nodes_prob <- nodes$values()[top_index]
      S <- sum(unlist(top_nodes_prob))
      top_nodes_prob <- lapply(top_nodes_prob, function(el) el / S)
      nodes <- dict(items = top_nodes_prob, keys = top_nodes)
    }
    
    if(store_mixture == TRUE) collection_nodes[[i]] <- nodes
    
    #partial.elapsed.time <- round(difftime(cur.time, aux.time), 3)
    #message(sprintf('Update in %.2f %s', partial.elapsed.time, units(partial.elapsed.time)))
    #elapsed.time <- round(difftime(sys.time(), start.time.PF), 3)
    #message('    Overall in ', elapsed.time, '', units(elapsed.time))
    
    #diagnostics: the distribution of the probabilities of active nodes
    #current_probs <- sort( unlist( nodes$values() ), decreasing = T )
    #plot(current_probs, type = 'h', main = 'Update', ylim = c(0, max(current_probs)) )
    #abline(h = pruning, lty = 'dashed', col = 'red')
    
    if (i <= length(observation_list)) {
      gamma <- observation_list[[i]]
      nodes <- UpdateUnnormalized2(nodes, gamma, alpha, theta, epsilon = epsilon, N = 0, show_top = Inf, pruning = pruning_update, verbose = verbose)[[1]]
    } 
    
    #diagnostics: the distribution of the probabilities of active nodes
    #current_probs <- sort( unlist( nodes$values() ), decreasing = T )
    #plot(current_probs, type = 'h', main = 'Propagation', ylim = c(0, max(current_probs)))
    #abline(h = pruning, lty = 'dashed', col = 'red')
    
  }
  
  elapsed.time <- round(difftime(Sys.time(), start.time.PF), 3)
  message('    Duality predictive runtime: ', elapsed.time, ' ', units(elapsed.time) )
  
  if(store_mixture == TRUE) return(list(predictive, collection_nodes)) else return(predictive)
}

#------------------------------------------------------------------------------
# Bootstrap particle filter
# Computes the likelihood (monomial symmetric function)
# WARNING! THE COMPUTATION OF THE LIKELIHOOD P_arrow() TAKES A LOT OF TIME! (apart from toy examples...)

P_arrow <- function(gamma, x){
  
  x <- x / sum(x)
  
  l_gamma <- length(gamma)
  
  if(l_gamma > length(x) ) return(0)
  
  x_subs <- permuteGeneral(x, m = l_gamma)
  
  normalization <- C(gamma)
  
  subs_sum <- apply(x_subs, 1, function(el) normalization*prod(el^gamma) )
  
  subs_sum <- sum(subs_sum)
  
  return(subs_sum)
  
}

Bootstrap_filter_King <- function(observation_list, alpha, theta, delta, epsilon, nParticles, show_top = 3){
  #DO NOT use this algorithm if delta is not constant
  #indeed this algorithm store and reuses the matrix B_km to speed up the computations
  
  print("-----------> PARTICLE BOOTSTRAP FILTER (Kingman's) <-----------")
  
  time.start <- Sys.time()
  
  filter <- vector( mode = 'list', length = length(observation_list) )
  importance_weights <- vector( mode = 'list', length = length(observation_list) )
  
  t <- 0
  
  update_size <- c()
  predict_size <- c()
  
  B_km <- NULL
  
  for( gamma in observation_list ){
    
    L <- length(gamma)
    t <- t + 1
    message(paste("Time", t, 'out of', length(observation_list)), appendLF = TRUE)
    
    #sample from stationary distribution if count == 1, else propagate
    if(t == 1) {
      message('Propagating...')
      empirical <- pblapply(1:nParticles, function(n) pd(alpha, theta, epsilon, show_top = Inf, min_length = L)$weights )
    }
    if(t > 1){
      #now resample (bootstrap)
      message('Resampling...')
      bootstrap <- sample(1:nParticles, nParticles, replace = TRUE, prob = weights)
      empirical <- pblapply( bootstrap, function(n) empirical[[n]] ) #this is an approximation of the filtering distribution
      
      #now propagate the resampled particles
      #message('Propagating...')
      #empirical <- pblapply( empirical, function(el) Petrov_OneStep(alpha, theta, delta = t, epsilon, X0 = el, show_top = Inf, min_length = L ) )
      temp <- Petrov_OneStep(alpha, theta, delta = delta, epsilon, X0 = empirical, show_top = Inf, min_length = L, B_km = B_km )
      empirical <- temp$empirical
      B_km <- temp$B_km
    }
    
    #now update: KINGMAN'S LIKELIHOOD
    message('Update...')
    weights <- unlist( pblapply( empirical, function(x) P_arrow(gamma, x)  ) )
    weights <- weights / sum(weights)
    

    #now save the first show_top for each atom
    #for (i in 1:show_top) {
      #top_i <- unlist( lapply(empirical, function(el) el[i] ) )
      #if( i == 1 ) top_t <- top_i else top_t <- cbind(top_t, top_i) #top_t is the matrix containing the top frequencies
    #}
    #filter[[count]] <- top_t
    
    #now save the first show_top for each atom
    #vector of frequencies + importance weights is an approximation of the filtering distribution!
    empirical <- lapply(empirical, function(x) {
      if (length(x) >= show_top) x[1:show_top] else c(x, rep(0, show_top - length(x)))
    })
    filter[[t]] <- do.call('rbind', empirical)
    importance_weights[[t]] <- weights

  }
  
  message('    Particle bootstrap filter (OLD) runtime (minutes): ', difftime(Sys.time(), time.start, units = 'mins') )
  return(list('Filter' = filter, 'Importance_weights' = importance_weights))
  
}

Bootstrap_filter_King_qentropy <- function(observation_list, alpha, theta, delta, epsilon, nParticles, MCrep_P_arrow = 10^4, show_top = Inf, q = 1){
  #WARNING! DO NOT use this algorithm if delta is not constant
  #indeed this algorithm store and reuses the matrix B_km to speed up the computations
  
  print("-----------> PARTICLE BOOTSTRAP FILTER (Kingman's) <-----------")
  
  time.start <- Sys.time()
  
  filter <- vector( mode = 'list', length = length(observation_list) )
  importance_weights <- vector( mode = 'list', length = length(observation_list) )
  
  t <- 0
  B_km <- NULL
  
  for( gamma in observation_list ){
    
    L <- length(gamma)
    t <- t + 1
    message(paste("Time", t, 'out of', length(observation_list)), appendLF = TRUE)
    
    #sample from stationary distribution if count == 1, else propagate
    if(t == 1) {
      message('    Propagating...')
      empirical <- pblapply(1:nParticles, function(n) pd(alpha, theta, epsilon, show_top = Inf, min_length = L)$weights )
    }
    if(t > 1) {
      #now resample (bootstrap)
      message('    Resampling...')
      bootstrap <- sample(1:nParticles, nParticles, replace = TRUE, prob = weights)
      empirical <- pblapply( bootstrap, function(n) empirical[[n]] ) #this is an approximation of the prediction distribution
      
      #message('    Propagating...')
      #empirical <- pblapply( empirical, function(el) Petrov_OneStep(alpha, theta, delta = t, epsilon, X0 = el, show_top = Inf, min_length = L ) )
      temp <- Petrov_OneStep(alpha, theta, delta = delta, epsilon, X0 = empirical, show_top = Inf, min_length = L, B_km = B_km )
      empirical <- temp$empirical
      B_km <- temp$B_km
    }
    
    #now update
    message('    Update...')
    weights <- pbsapply( empirical, function(x) P_arrow(gamma, x) ) #without parallelization
    #weights <- pbmclapply( empirical, function(x) P_arrow(gamma, x), mc.cores = detectCores() - 1 ) #with parallelization
    #weights <- unlist(weights) #with parallelization
    weights <- weights / sum(weights)
    
    #now save the zigosity/entropy of order q for each atom
    #vector of frequencies + importance weights is an approximation of the filtering distribution!
    filter[[t]] <- sapply(empirical, function(x) qentropy(x / sum(x), q) )
    importance_weights[[t]] <- weights
    
  }
  
  message('    Particle bootstrap filter runtime (minutes): ', difftime(Sys.time(), time.start, units = 'mins') )
  
  return(list('Filter' = filter, 'Importance_weights' = importance_weights))
  
}

#-------------------------------------------------------------------------------
# Smoothing and backward-information recursion

Smoothing_one <- function(observation_list, target, alpha, theta, t, epsilon, pruning = 0, N = 10^4, nParticles = 10^4, show_top = 3, verbose = FALSE, bootstrap = FALSE){
  
  steps <- length(observation_list)
  forward_obs <- lapply(1:target, function(i) observation_list[[i]] )
  backward_obs <- lapply( steps:(target+1), function(i) observation_list[[i]] )
  print(backward_obs)
  
  print('-----------> FORWARD FILTER <-----------')
  filter <- Particle_filter( forward_obs, alpha, theta, t, epsilon, pruning, N, nParticles, show_top, verbose, bootstrap )
  
  forward_nodes <- filter[[4]]
  filter <- filter[[1]]
  
  print('-----------> COST TO GO <-----------')
  backward_nodes <- Cost_to_go( backward_obs , alpha, theta, t, epsilon, pruning, nParticles, verbose, bootstrap )
  backward_nodes <- backward_nodes[[1]]
  
  print('-----------> SMOOTHING <-----------')
  random_obs <- runif( N )
  nodes_distribution <- cumsum( unlist( backward_nodes$values() ) / nParticles )
  nodes_index <- unlist( lapply( random_obs, function(x) which.max( x <= nodes_distribution ) ) )
  random_obs <- lapply( nodes_index, function(i) backward_nodes$keys()[[i]] )
  
  random_start <- runif( N )
  nodes_distribution <- cumsum( unlist( forward_nodes$values() ) / nParticles )
  nodes_index <- unlist( lapply( random_start, function(x) which.max( x <= nodes_distribution ) ) )
  random_start <- lapply( nodes_index, function(i) forward_nodes$keys()[[i]] )
  
  for( i in 1:N ){
    
    status <- i / N * 100
    if( status %% 5 == 0 ) print(paste('    ...status:', round( status, 2), '%' ) )
    
    current_update <- CoagulationSet( random_start[[i]], random_obs[[i]], alpha, theta, verbose = F )
    current_prob <- lapply( current_update$values() , function(el) el$prob )
    current_sample <- sample( 1:current_update$size() , 1, prob = current_prob  )
    current_node <- current_update$keys()[[current_sample]]
    
    random_measure <- pd_posterior(alpha, theta, current_node, epsilon)$weights[1:show_top]
    
    if(i == 1) smoothing <- random_measure else smoothing <- rbind(smoothing, random_measure)
    
  }
  
  print('    ...done!')
  
  return( list( filter, smoothing ) )
  
}

Cost_to_go <- function( observation_list, alpha, theta, t, epsilon, pruning_update = 0, pruning_propagation = 0, nParticles = 10^4, verbose = FALSE ){
  
  #N is the number of Monte Carlo replicates per observed partition
  start.time.PF <- Sys.time()
  L <- length(observation_list)
  nodes <- dict( keys = list(observation_list[[L]]), items = nParticles)
  collection_nodes <- vector(mode = 'list', length = length(observation_list) )
  
  count <- 0
  
  for(i in (L-1):1 ){
    
    #aux.time <- Sys.time()
    
    count <- count + 1
    message(paste("From", i + 1, 'to', i ) )
    
    nodes <- Particle_predict(nodes, t, theta, pruning = pruning_propagation, nParticles = nParticles, verbose = verbose )
    collection_nodes[[i]] <- nodes
    
    #partial.elapsed.time <- round(difftime(cur.time, aux.time),3)
    #message(sprintf('Propagation in %.2f %s', partial.elapsed.time, units(partial.elapsed.time)))
    #elapsed.time <- round(difftime(Sys.time(), start.time.PF),3)
    #message('    Runtime', elapsed.time, '', units(elapsed.time))
    
    if(length(nodes$keys()) == 0){
      print('Attenzione! Pruning troppo alto! Tutti i nodi sono stati eliminati!')
      return(NULL)
    }
    
    #diagnostics: the distributions of the probabilities of active nodes
    #current_probs <- sort( unlist( nodes$values() ), decreasing = T )
    #plot(current_probs, type = 'h', main = 'Update', ylim = c(0, max(current_probs)) )
    #abline(h = pruning, lty = 'dashed', col = 'red')
    
    #aux.time <- Sys.time()
    
    gamma <- observation_list[[i]]
    if(i > 1) nodes <- UpdateUnnormalized2(nodes, gamma, alpha, theta, N = 0, pruning = pruning_update, verbose = verbose)[[1]]
    
    #partial.elapsed.time <- round(difftime(cur.time,aux.time),3)
    #message(sprintf('Update in %.2f %s', partial.elapsed.time, units(partial.elapsed.time)))
    #elapsed.time <- round(difftime(Sys.time(), start.time.PF),3)
    #if(i > 1) message( '    Overall in %.2f %s', elapsed.time, ' ', units(elapsed.time))
    
    # print(paste('Active nodes post pruning:', nodes$size()))
    
    #diagnostics: the distribution of the probabilities of active nodes
    #current_probs <- sort( unlist( nodes$values() ), decreasing = T )
    #plot(current_probs, type = 'h', main = 'Propagation', ylim = c(0, max(current_probs)))
    #abline(h = pruning, lty = 'dashed', col = 'red')
    
  }
  
  #plot(update_size, type = 'b', pch = 16, main = paste('Update set size. Pruning=', pruning) )
  #plot(predict_size, type = 'b', pch = 16, main = paste('Propagation set size. Pruning=', pruning) )
  
  elapsed.time <- round(difftime(Sys.time(), start.time.PF), 3)
  message('    Backward filter runtime: ', elapsed.time, ' ', units(elapsed.time))
  
  return( collection_nodes )
}

Smoothing <- function( observation_list, alpha, theta, t, epsilon, pruning_update = 0, pruning_propagation = 0, N = 10^4, nParticles = 10^4, show_top = 3, verbose = FALSE, bootstrap = FALSE){
  
  time.start <- Sys.time()
  steps <- length(observation_list)
  
  message('-----------> FORWARD FILTER <-----------')
  filter <- Particle_filter( observation_list, alpha, theta, t, epsilon, pruning_update = pruning_update, pruning_propagation = pruning_propagation, N, nParticles, show_top, verbose, store_mixture_update = TRUE )
  forward_nodes <- filter[[2]]
  filter <- filter[[1]]
  
  message('-----------> BACKWARD FILTER <-----------')
  backward_nodes <- Cost_to_go( observation_list, alpha, theta, t, epsilon, pruning_update = pruning_update, pruning_propagation = pruning_propagation, nParticles, verbose )
  
  message('-----------> SMOOTHING <-----------')
  smoothing <- vector(mode = 'list', length = steps-1 )
  
  if ( (pruning_update < 1 && pruning_update > 0.1) ) { #|| pruning_update > 10 ) {
    for(target in 1:(steps-1) ){
    
    message("Time ", target, ' out of ', steps - 1 )
    
    current_forward <- forward_nodes[[target]]
    forward_probs <- unlist( current_forward$values() )
    index <- sample(1:current_forward$size(), N, replace = TRUE, prob = forward_probs )
    #random_start <- lapply( index, function(i) current_forward$keys()[[i]] )
    random_start <- current_forward$keys()[index]
    
    current_backward <- backward_nodes[[ target ]]
    backward_probs <- unlist( current_backward$values() )
    psf_backward <- sapply(current_backward$keys(), function(el) PSF(alpha, theta, el, normalized = TRUE) )
    backward_probs <- backward_probs / psf_backward
    index <- sample(1:current_backward$size(), N, replace = TRUE, prob = backward_probs)
    #proposals <- lapply(index, function(i) current_backward$keys()[[i]])
    proposals <- current_backward$keys()[index]
    
    message('    Coagulating forward and backward filters (urns)...')
    mus <- pblapply(1:N, function(i) Polya_urn(alpha, theta, random_start[[i]], proposals[[i]] ) )
    smoothing[[target]] <- t(sapply(mus, function(mu) pd_posterior(alpha, theta, mu, epsilon, show_top = show_top)$weights ))
    
    }
  }
  else {
    for(target in 1:(steps-1) ){
    
    message("Time ", target, ' out of ', steps - 1 )
    
    message('    Coagulating forward and backward filters (exact)...')
    current_forward <- forward_nodes[[target]]
    current_backward <- backward_nodes[[target]]
    
    coag_set <- dict()
    
    pblapply(current_backward$keys(), function(gamma)
      lapply(current_forward$keys(), function(omega) {
        temp <- CoagulationSetUnnormalized(omega, gamma, alpha, theta, verbose = F)
        prob_omega <- current_forward$get(omega)
        prob_gamma <- current_backward$get(gamma)
        
        for(mu in temp$keys() ){
          
          prob_old <- coag_set$get(as.numeric(mu), default = 0 )
          prob_new <- prob_old + temp$get(mu)$prob * prob_omega / PSF(alpha, theta, omega, normalized = FALSE) * prob_gamma / PSF(alpha, theta, gamma, normalized = FALSE)
          
          coag_set$set(as.numeric(mu), prob_new )
          
        }
      
      })
    )
    
    S <- sum(unlist(coag_set$values() ))
    coag_set <- dict(keys = coag_set$keys(), items = as.list(unlist(coag_set$values()) / S) )
    
    message('    Drawing from smoothing distribution...')
    current_sample <- sample(coag_set$keys(), N, replace = TRUE, prob = unlist(coag_set$values()) )
    smoothing[[target]] <- t(pbsapply(current_sample, function(mu) pd_posterior(alpha, theta, mu, epsilon, show_top = show_top)$weights ))
    }
    
    smoothing[[steps]] <- filter[[steps]]
  }
  
  
return( list(filter, smoothing) )
  
}

Smoothing_qentropy <- function( observation_list, alpha, theta, t, epsilon, pruning_update = 0, pruning_propagation = 0, N = 10^4, nParticles = 10^4, show_top = 3, verbose = FALSE, bootstrap = FALSE, q = 1){
  
  time.start <- Sys.time()
  steps <- length(observation_list)
  
  message('-----------> FORWARD FILTER <-----------')
  filter <- Particle_filter_qentropy( observation_list, alpha, theta, t, epsilon, pruning_update = pruning_update, pruning_propagation = pruning_propagation, N, nParticles, show_top, verbose, store_mixture = TRUE, q = q )
  forward_nodes <- filter[[2]]
  filter <- filter[[1]]
  
  message('-----------> BACKWARD FILTER <-----------')
  backward_nodes <- Cost_to_go( observation_list, alpha, theta, t, epsilon, pruning_update = pruning_update, pruning_propagation = pruning_propagation, nParticles, verbose )
  
  message('-----------> SMOOTHING <-----------')
  smoothing <- vector(mode = 'list', length = steps )
  
  if ( (pruning_update < 1 && pruning_update > 0.1) ) { #|| pruning_update > 10 ) {
    for(target in 1:(steps-1) ){
      
      message("Time ", target, ' out of ', steps - 1 )
      
      current_forward <- forward_nodes[[target]]
      forward_probs <- unlist( current_forward$values() )
      index <- sample(1:current_forward$size(), N, replace = TRUE, prob = forward_probs )
      #random_start <- lapply( index, function(i) current_forward$keys()[[i]] )
      random_start <- current_forward$keys()[index]
      
      current_backward <- backward_nodes[[ target ]]
      backward_probs <- unlist( current_backward$values() )
      psf_backward <- sapply(current_backward$keys(), function(el) PSF(alpha, theta, el, normalized = TRUE) )
      backward_probs <- backward_probs / psf_backward
      index <- sample(1:current_backward$size(), N, replace = TRUE, prob = backward_probs)
      #proposals <- lapply(index, function(i) current_backward$keys()[[i]])
      proposals <- current_backward$keys()[index]
      
      message('    Coagulating forward and backward filters (urns)...')
      mus <- pblapply(1:N, function(i) Polya_urn(alpha, theta, random_start[[i]], proposals[[i]] ) )
      smoothing[[target]] <- t(sapply(mus, function(mu) pd_posterior(alpha, theta, mu, epsilon, show_top = show_top)$weights ))
      
    }
  }
  else {
    for(target in 1:(steps-1) ){
      
      message("Time ", target, ' out of ', steps - 1 )
      
      message('    Coagulating forward and backward filters (exact)...')
      current_forward <- forward_nodes[[target]]
      current_backward <- backward_nodes[[target]]
      
      coag_set <- dict()
      
      lapply(current_backward$keys(), function(gamma) #use pblapply
        lapply(current_forward$keys(), function(omega) {
          temp <- CoagulationSetUnnormalized(omega, gamma, alpha, theta, verbose = F)
          prob_omega <- current_forward$get(omega)
          prob_gamma <- current_backward$get(gamma)
          
          for(mu in temp$keys() ){
            
            prob_old <- coag_set$get(as.numeric(mu), default = 0 )
            prob_new <- prob_old + temp$get(mu)$prob * prob_omega / PSF(alpha, theta, omega, normalized = FALSE) * prob_gamma / PSF(alpha, theta, gamma, normalized = FALSE)
            
            coag_set$set(as.numeric(mu), prob_new )
            
          }
          
        })
      )
      
      S <- sum(unlist(coag_set$values() ))
      coag_set <- dict(keys = coag_set$keys(), items = as.list(unlist(coag_set$values()) / S) )
      
      message('    Drawing from smoothing distribution...')
      current_sample <- sample(coag_set$keys(), N, replace = TRUE, prob = unlist(coag_set$values()) )
      temp <- lapply(current_sample, function(mu) pd_posterior(alpha, theta, mu, epsilon, show_top = show_top)$weights ) #use pblapply
      smoothing[[target]] <- sapply(temp, function(el) qentropy(el / sum(el), q = q) )
    }
    
    smoothing[[steps]] <- filter[[steps]]
  }
  
  message('    Dual smoother runtime (minutes): ', difftime(Sys.time(), time.start, units = 'mins') )
  
  return( list(filter, smoothing) )
  
}

