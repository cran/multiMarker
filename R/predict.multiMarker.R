
#--------------#
#- Prediction -#
#--------------#


predict.multiMarker <- function( object, y,
                                 niter = 10000, burnIn = 3000,
                                 posteriors = FALSE, ...){

  burn <- seq( burnIn, niter)
  yNew <- y
  P <- object$constants$P
  tmp_P <- dim( yNew)
  if( is.null(tmp_P )){
    tmp_P <- 1
    checkP <- length(yNew)
  }else{
    nNew <- tmp_P[1]
    tmp_P <- tmp_P[2]
    checkP <- tmp_P
  }
  # one i at the time
  try( if( checkP != P ) stop(
    "Number of metabolites in test data different from number of metabolites in train data!"))

  #--- initialize ---#
  D <- object$constants$D
  alphaUp <- object$estimates$ALPHA_E[1,]
  alphaUpsd <- object$estimates$ALPHA_E[2,]
  alphaUpT <- alphaUp
  betaUp <- object$estimates$BETA_E[1,]
  betaUpsd <- object$estimates$BETA_E[2,]
  betaUpT <- betaUp
  sigma2_errUp <- object$estimates$SigmaErr_E[1,]
  varDUp <- object$estimates$SigmaD_E[1,]
  x_D <- unique(object$constants$x_D)
  tauD <- object$constants$tauD
  thetaUp <- object$estimates$THETA_Est[,,1]
  thetasd <-  object$estimates$THETA_Est[,,2]
  boundsL <- c(-Inf, thetaUp[1,-(D-1)])
  boundsU <- c(thetaUp[1,-1], Inf)
  s2EHp <- object$estimates$varPHp

  #-- old hps --#
  thetaM <-  thetaUp  # (P+1) x D-1
  sigmaWprior <- object$constants$sigmaWprior

  if( tmp_P == 1){
    #--- probabilities ---#
    probs_c <- cauchit_probs(yNew, thetaUp, D) # vector of length D
    labelNew <- which.max(probs_c)
    #--- intake initialization ---#
    zPred <- x_D[labelNew] + rnorm(1, 0, 1)
    tmpNA <- which(is.na(yNew) )
    if( length(tmpNA) != 0 ){
      for(kk in tmpNA ){
        yNew[kk] <- object$constants$y_Median[kk]
      }
    }
    #-- STORING --#
    ZPRED <- rep( NA, niter)
    ZPRED[1] <- zPred
    PROBS <- matrix(NA, nrow = niter, ncol = D)
    PROBS[1, ] <- probs_c
  }else{
    probs_c <- t(apply(yNew, 1, function(x) cauchit_probs(x, thetaUp, D) ))# matrix with n rows and D cols
    labelNew <- apply(probs_c, 1, which.max)
    #--- intake initialization ---#
    zPred <- x_D[labelNew] + rnorm(nNew, 0, 1)
    tmpNA <- which(is.na(yNew) )
    if( length(tmpNA) != 0 ){
      for(kk in tmpNA ){
        where_temp <- floor((kk-1)/nNew) +1
        yNew[kk] <- object$y_Median[where_temp]
      }
    }
    #-- STORING --#
    ZPRED <- matrix( NA, ncol = niter, nrow = nNew)
    ZPRED[,1] <- zPred
    PROBS <- array(NA, dim = c(nNew, D, niter ))
    PROBS[,,1] <- probs_c
  }

  pbar <- txtProgressBar(min = 2, max = (niter + 1), style = 3)
  on.exit(close(pbar))

  #--- iterations ---#

  for ( it in 2:niter){

    setTxtProgressBar(pbar, it)

    # alpha
    alphaUp <- sapply(1:P, function(p)
      rtruncnorm(1, 0, Inf,  alphaUpT[p], alphaUpsd[p]) ) #

    # # beta
    betaUp <- sapply(1:P, function(p)
      rtruncnorm(1, 0, Inf, betaUpT[p], betaUpsd[p] ) )  #
    #
    # # sigma2P
    sigma2_errUp <- sapply(1:P, function(p)
      1/rgamma(1, shape = s2EHp[1,p], rate = s2EHp[2,p])  )

    # # theta
    thetaUp[1,] <- sapply(1:(D-1), function(d)
      rtruncnorm(1, boundsL[d], boundsU[d], thetaM[1,d], thetasd[1,d] ) )

    thetaUp[-1,] <- apply( thetaM[-1,], c(1,2), function(x)
      rnorm( 1, x, 0.0001))

    # compute new probs and sample new z value/values
    if( tmp_P == 1){

      probs_c <- cauchit_probs(yNew, thetaUp, D)
      labelNew <- which.max( probs_c)
      zPred <- z_fc(  varDUp[labelNew], x_D[labelNew], sigma2_errUp,
                      betaUp, yNew, alphaUp, P, tauD[labelNew] )
      ZPRED[it] <- zPred
      PROBS[it, ] <- probs_c
    }else{
      probs_c <- t(apply(yNew, 1, function(x) cauchit_probs(x, thetaUp, D) ))# matrix with n rows and D cols
      labelNew <- apply(probs_c, 1, which.max)
      zPred <- sapply( 1:nNew, function(i)
        z_fc(  varDUp[labelNew[i]], x_D[labelNew[i]], sigma2_errUp,
               betaUp, yNew[i,], alphaUp, P, tauD[labelNew[i]] ))
      ZPRED[,it] <- zPred
      PROBS[,,it] <- probs_c
    }
  } # close niter loop

  if( tmp_P == 1){
    ZPred_P <- msdci(ZPRED[burn])
    Prob_P1 <- apply(PROBS[burn,], 2, median)
    Prob_P2 <- apply(PROBS[burn,], 2, sd)
    Prob_P3 <- apply(PROBS[burn,], 2, function(x) quantile(x, 0.025))
    Prob_P4 <- apply(PROBS[burn,], 2, function(x) quantile(x, 0.975))
    Prob_P <- rbind( Prob_P1,  Prob_P2,  Prob_P3,  Prob_P4)

    predictions <- list( ZPred_P = ZPred_P, Prob_P = Prob_P)
  }else{
    ZPred_P <- sapply(1:nNew, function(i) msdci(ZPRED[i,burn]))
    Prob_P1 <- apply(PROBS[,,burn], c(1,2), median)
    Prob_P2 <- apply(PROBS[,,burn], c(1,2), sd)
    Prob_P3 <- apply(PROBS[,,burn], c(1,2), function(x) quantile(x, 0.025))
    Prob_P4 <- apply(PROBS[,,burn], c(1,2), function(x) quantile(x, 0.975))
    Prob_P <- array( cbind( Prob_P1,  Prob_P2,  Prob_P3,  Prob_P4), dim = c(nNew, D, 4))

    predictions <- list( ZPred_P = ZPred_P, Prob_P = Prob_P)
  }

  chains <- list( ZPRED = ZPRED, PROBS = PROBS)

  out <- list( predictions = predictions,
               chains = if(posteriors) chains else NULL)
  return(out)

} # close function
