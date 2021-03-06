##################################################
## multiMarker model - Estimation + PRrediction ##
##################################################


multiMarker <- function(y, quantities,
                        niter = 10000, burnIn = 3000,
                        posteriors = FALSE, sigmaAlpha = 1,
                        nuZ1 = NULL, nuZ2 = NULL,
                        nuSigmaP1 = NULL, nuSigmaP2 = NULL, sigmaWprior = 0.000001,
                        nuBeta1 = 2, nuBeta2 = 3, tauBeta = 0.1
){

  call <- match.call()

  #--- constant quantities ---#
  n <- nrow(y)
  P <- ncol(y)
  temp <- table(quantities)
  D <- dim(temp)
  labels_D <- as.numeric(as.factor(quantities))
  x_D <- sort(as.numeric( names(temp)))
  burn <- seq( burnIn, niter)

  n_D <- as.numeric(temp) # number of subjects per component
  pi_D <- n_D/sum(n_D)

  #--- initialization ---#
  tmp_lm <- lm(y ~ quantities)
  mAlpha  <- mean( tmp_lm$coefficients[1,])
  mBeta <-  mean( tmp_lm$coefficients[2,])

  tempAB <- alpha_beta_iniz( x_D, y, D, n_D, P)
  alpha <-  tempAB$alpha_iniz
  beta <-  tempAB$beta_iniz
  tauAlpha = 1


  tempS <- sigma2_err_iniz(beta, y, D, P)
  sigma2_err <- sapply(1:P, function(x) summary(tmp_lm)[[x]]$sigma *(1/P))^2
  if( is.null(nuSigmaP1) ){nuSigmaP1 <- 3}
  if( is.null(nuSigmaP2) ){nuSigmaP2 <- 1}
  if( is.null(nuZ2) ){
    nuZ2 <- rep(n, D)
  }
  if( is.null(nuZ1) ){
    tmp <- diff(c(0, x_D))
    nuZ1 <- seq(D,1, length = D)/2
    nuZ1[ tmp > 2*median(tmp) ] <- nuZ1[ tmp > 2*median(tmp) ]*2
  }

  muAlpha <- mean(alpha); muBeta <- mean(beta); varBeta <- var(beta)
  varD <- rep(5^2, D)


  z <- rtruncnorm(n, 0, Inf, quantities, sqrt(varD))
  tmp <- cbind(labels_D, z)
  muD <- x_D

  leva <- which(is.na(rowSums(y)))
  if( length(leva) == 0 ){
    mod0 <- ordinalNet(as.matrix(y) , as.factor(quantities),
                       family = "cumulative", link ="cauchit" ,
                       parallelTerms=TRUE)
  }else{
    mod0 <- ordinalNet(as.matrix(y)[-leva,] , as.factor(quantities)[-leva],
                       family = "cumulative", link ="cauchit" ,
                       parallelTerms=TRUE)
  }
  theta <- coef(mod0, matrix = T)
  thetaM <- theta
  probs_c <- t( apply(y, 1, function(x) cauchit_probs(x, theta, D)) )
  PROBS <- array(NA, dim = c(n, D, niter ))
  PROBS[,,1] <- probs_c

  #--- storing ---#
  ALPHA <- BETA <- SigERR <- matrix(NA, nrow = niter, ncol = P)
  Z <- matrix(NA, nrow = niter, ncol = n)
  ALPHABETAPAR <- matrix(NA, nrow = niter, ncol = 3) # (in order) sigmaBeta, muAlpha, muBeta
  MUD <- SIGD <- matrix(NA, nrow = niter, ncol = D)
  STORE_hp_s2P <- matrix(NA, nrow = 2, ncol = P)
  tracker <- 1
  ALPHA[1, ] <- alpha; BETA[1, ] <- beta; SigERR[1, ] <- sigma2_err
  Z[1, ] <- z
  ALPHABETAPAR[1, ] <- c(varBeta, muAlpha, muBeta)
  MUD[1, ] <- muD; SIGD[1, ] <- varD
  THETA <- array(NA, dim = c(P+1, D -1, niter)) # first row is class- intercepts, columns are P scaling parameters (constant in d=1,..,D-1)
  ACCTHETA <- array(NA, dim = c(P+1, D-1, niter))
  THETA[,,1] <- theta
  thetaM <- theta
  boundsL <- c(-Inf, theta[1,-(D-1)])
  boundsU <- c(theta[1,-1], Inf)
  y_Median <- apply(y, 2, function(x) median(x, na.rm = T))
  y_Var <- apply(y, 2, function(x) var(x, na.rm = T))

  #--- compute scale factor ---#

  labelNew <- t(sapply(1:n,
                       function(i) sapply( 1:D,
                                           function(d)
                                             probs_c[i,d]*dtruncnorm(z[i], 0, Inf, x_D[d], sqrt(varD[d])) )))
  labelNew <- labelNew/rowSums(labelNew)
  labelNew <- apply(labelNew,1,which.max)


  #------#

  pbar <- txtProgressBar(min = 2, max = (niter + 1), style = 3)
  on.exit(close(pbar))

  #--- MCMC ---#

  for (it in 2:niter){

    setTxtProgressBar(pbar, it)

    # UPDATE alpha-beta prior parameters
    varBeta <- variance_fc(beta, P, nuBeta1, nuBeta2,
                           tauBeta, muBeta, mBeta)
    muAlpha <- mean_fc( tauAlpha, sigmaAlpha, alpha, P,
                        mAlpha)
    muBeta <- mean_fc( tauBeta, varBeta, beta, P,
                       mBeta)
    ALPHABETAPAR[it, ] <- c(varBeta, muAlpha, muBeta)

    # UPDATE alpha-beta
    alphaP <- sapply(1:P, function(p) alpha_fc( sigmaAlpha, beta[p], n, muAlpha,
                                                z, sigma2_err[p], y[,p]) )
    alpha <- unlist(alphaP[1,])
    ALPHA[it,] <- alpha

    betaP <- sapply(1:P, function(p) beta_fc( varBeta, alpha[p], n, muBeta,
                                              z, sigma2_err[p], y[,p]))
    beta <- unlist(betaP[1,])
    BETA[it,] <- beta

    # update error variances
    sigma2_errP <- sapply(1:P, function(p) sigma2_err_fc( n, nuSigmaP1, nuSigmaP2,
                                                          y[,p],
                                                          alpha[p], beta[p], z))
    sigma2_err <- unlist(sigma2_errP[1,])

    SigERR[it, ] <- sigma2_err

    # update component parameters
    n_D <- tabulate( as.factor(sort(labelNew)))

    varDP <- sapply(1:D, function(d) variance_fc_d( z[which(labelNew == d)],
                                                    n_D[d], nuZ1[d], nuZ2[d],
                                                    1, x_D[d],
                                                    0))

    varD <- unlist(varDP[1,])

    SIGD[it, ] <- varD

    # update intercept param for weights params
    int_prop <- rep(NA, D-1)
    bl_prop <- boundsL
    bu_prop <- boundsU
    for ( d in 1:(D-1)){
      int_prop[d] <- rtruncnorm(1, boundsL[d], boundsU[d], thetaM[1,d], sigmaWprior )
    }

    thetaTemp <- theta
    thetaTemp[1,] <- int_prop
    probs_cTemp <- t( apply(y, 1, function(x) cauchit_probs(x, thetaTemp, D)) )
    lAccTheta <- logPost_MCMC_wcum(probs_cTemp, D, n, labelNew) -
      logPost_MCMC_wcum(probs_c, D, n, labelNew)

    logU <- log(runif(1))
    indAcc <- ( lAccTheta < logU )
    indAcc[is.nan(lAccTheta)] <- TRUE; indAcc[is.infinite(lAccTheta)] <- TRUE
    if( is.nan(indAcc)){indAcc <- TRUE}

    if(  ( indAcc  ) ){ # reject
      ACCTHETA[1,, it] <- 0
    } else { # accept
      ACCTHETA[1,, it] <- 1
      theta <- thetaTemp # update values
      probs_c <- probs_cTemp
    }

    # update coeff param for weights params
    int_prop <- apply( thetaM[-1,], c(1,2), function(x)
      rnorm( 1, x, sigmaWprior))
    thetaTemp <- theta
    thetaTemp[-1,] <- int_prop
    probs_cTemp <- t( apply(y, 1, function(x) cauchit_probs(x, thetaTemp, D)) )
    lAccTheta <- logPost_MCMC_wcum(probs_cTemp, D,n, labelNew) -
      logPost_MCMC_wcum(probs_c, D, n, labelNew)

    logU <- log(runif(1))
    indAcc <- ( lAccTheta < logU )
    indAcc[is.nan(lAccTheta)] <- TRUE; indAcc[is.infinite(lAccTheta)] <- TRUE
    if( is.nan(indAcc)){indAcc <- TRUE}

    if(  ( indAcc  )  ){ # reject
      ACCTHETA[-1,, it] <- 0
    } else { # accept
      ACCTHETA[-1,, it] <- 1
      theta <- thetaTemp # update values
      probs_c <- probs_cTemp
    }
    THETA[,,it] <- theta

    # compute new probs and update latent intakes

    labelNew <- t(sapply(1:n,
                         function(i) sapply( 1:D,
                                             function(d)
                                               probs_c[i,d]*dtruncnorm(z[i], 0, Inf, x_D[d], sqrt(varD[d])) )))
    labelNew2 <- labelNew/rowSums(labelNew)
    labelNew <- apply(labelNew2 , 1, function(x) sample(seq(1,D), 1, prob = x) )
    labelNew[which(is.na(labelNew))] <- sample(seq(1,D), length(which(is.na(labelNew))), replace = TRUE)

    z <- sapply( 1:n, function(i)
      z_fc(  varD[labelNew[i]], x_D[labelNew[i]], sigma2_err,
             beta, y[i,], alpha, P, 1 ))

    Z[it, ] <- z
    PROBS[,,it] <- probs_c

    #---#
    if( it == burnIn[1]){
      STORE_hp_s2P[1,] <- unlist(sigma2_errP[2,])
      STORE_hp_s2P[2,] <- unlist(sigma2_errP[3,])
      tracker <- tracker +1
    }
    if( it > burnIn[1]){
      STORE_hp_s2P[1,] <- (tracker *  STORE_hp_s2P[1,] +unlist(sigma2_errP[2,]))/(tracker +1)
      STORE_hp_s2P[2,] <- (tracker *  STORE_hp_s2P[2,] +unlist(sigma2_errP[3,]))/(tracker +1)
      tracker <- tracker +1
    }
  }

  #--- OUTPUT ---#
  #-- chains
  sigmaBeta_c <- ALPHABETAPAR[burn,1]
  muAlpha_c <- ALPHABETAPAR[burn,2]
  muBeta_c <- ALPHABETAPAR[burn,3]
  ALPHA_c <- ALPHA[burn, ]
  BETA_c <- BETA[burn, ]
  SigmaErr_c <- SigERR[burn, ]
  SigmaD_c <- SIGD[burn, ]
  Z_c <- Z[burn, ]
  THETA_c <- THETA[,,burn]

  #-- estimates
  sigmaBeta_E <- msdci(sigmaBeta_c)
  muAlpha_E <- msdci(muAlpha_c)
  muBeta_E <- msdci(muBeta_c)
  ALPHA_E <- sapply(1:P, function(p) msdci(ALPHA_c[,p]))
  BETA_E <- sapply(1:P, function(p) msdci(BETA_c[,p]))
  SigmaErr_E <- sapply(1:P, function(p) msdci(SigmaErr_c[,p]))
  SigmaD_E <- sapply(1:D, function(d) msdci(SigmaD_c[,d]))
  Z_E <- sapply(1:n, function(i) msdci(Z_c[,i]))

  THETA_E <- array(NA, dim =c(nrow(THETA_c), ncol(THETA_c),4))
  THETA_E[,,1] <- apply(THETA_c, c(1,2), median)
  THETA_E[,,2] <- apply(THETA_c, c(1,2), sd)
  THETA_E[,,3] <- apply(THETA_c, c(1,2), function(x) quantile(x,0.25))
  THETA_E[,,4] <- apply(THETA_c, c(1,2), function(x) quantile(x,0.75))

  #-- extra information
  acc_probs <- apply(ACCTHETA[,,-1], c(1,2), sum)/dim(ACCTHETA)[3]
  weights_info <- list(acc_probs = acc_probs)

  #-- output
  constants <- list(nuZ1 = nuZ1, nuZ2 = nuZ2,
                    sigmaAlpha = sigmaAlpha,
                    nuSigmaP1 = nuSigmaP1, nuSigmaP2 = nuSigmaP2,
                    nuBeta1 = nuBeta1, nuBeta2 = nuBeta2,
                    tauBeta = tauBeta, x_D = x_D, P = P, D = D, n = n,
                    sigmaWprior = sigmaWprior, y_Median = y_Median, y_Var = y_Var)

  chains <- list( ALPHA_c = ALPHA_c, BETA_c = BETA_c, SigmaErr_c = SigmaErr_c,
                  SigmaD_c = SigmaD_c, Z_c = Z_c, THETA_c = THETA_c,
                  sigmaBeta_c = sigmaBeta_c, muAlpha_c = muAlpha_c,
                  muBeta_c = muBeta_c, weights_info = weights_info)

  estimates <- list( ALPHA_E = ALPHA_E, BETA_E = BETA_E, SigmaErr_E = SigmaErr_E,
                     SigmaD_E = SigmaD_E, Z_E = Z_E, THETA_Est = THETA_E,
                     sigmaBeta_E = sigmaBeta_E, muAlpha_E = muAlpha_E,
                     muBeta_E = muBeta_E, varPHp = STORE_hp_s2P
  )

  out <- list(estimates = estimates, constants = constants,
              chains = if(posteriors) chains else NULL)

  class(out) <- "multiMarker"

  return(out)

}


print.multiMarker <- function(x, ...)
{
  P <- x$constants$P
  D <- x$constants$D

  h1 <- paste("multiMarker model with", P, "biomarkers and", D, "portions", sep = " ")
  sep <- paste0( rep("=", 27 + floor(P/10) + floor(D/10), collapse = "" ))
  cat("\n", sep, "\n")
  cat("  ", h1, "\n")
  cat("", sep, "\n", "\n")
}

