##--- utilities ---##
msdci <- function(p_chain){ # after burn in

  out <- c( median(p_chain), sd(p_chain),
            as.numeric( quantile(p_chain, c(0.025, 0.975), na.rm = TRUE) ) )
  return(out)
}


##############################
## initialization functions ##
##############################

alpha_beta_iniz <- function( x_D, y, D, n_D, P){

  tmpN <- c(1,cumsum(n_D))
  meanY_PD <-  sapply(1:P, function(p) # P cols and D rows
    sapply(1:D, function(d)
      mean(y[tmpN[d]:tmpN[d+1], p], na.rm = TRUE)
    ))

  ncolM <- (D-1)*D - sum( seq(1,D-1) )
  A_Mat <- B_Mat <- matrix( nrow = P, ncol = ncolM )

  # computing alpha and beta values
  for( p in 1:P){
    a <- 1
    for( d in D:2){
      for( j in (d-1):1){
        alpha <- max(0, ( meanY_PD[d,p] - meanY_PD[j,p]*(x_D[d]/x_D[j]) )/( 1 - (x_D[d]/x_D[j])) )
        beta <- (meanY_PD[j,p] - alpha)/x_D[j]
        if( (beta <= 0) | is.nan(beta)){beta <- 0.001}
        A_Mat[p, a] <- alpha
        B_Mat[p, a] <- beta
        a <- a+1
      }
    }
  } # close P loop (rows are the Ps)

  alpha_iniz <- rowMeans(A_Mat, na.rm = TRUE)
  beta_iniz <- rowMeans(B_Mat, na.rm = TRUE)

  out <- list( alpha_iniz =  alpha_iniz,  beta_iniz =  beta_iniz)
  return(out)
}


sigma2_err_iniz <- function(beta, y, D, P){

  temp1 <- diag( beta%*%t(beta) )/D
  temp2 <- diag( var(y) )

  sigma2_err_iniz <- rep( sum( temp2 - temp1, na.rm = TRUE ), P )
  sigma2_err_iniz[which(sigma2_err_iniz <= 0)] <- 0.3

  out <- list( sigma2_err_iniz  = sigma2_err_iniz )
  return(out)
}

z_par_iniz <- function(y, sigma2_err, sigmaAlpha,
                       sigmaBeta, P, D, n_D){
  tmpN <- c(1,cumsum(n_D))

  temp <- matrix(NA, ncol = P, nrow = D)
  varY_PD <-  sapply(1:P, function(p) # P cols and D rows
    sapply(1:D, function(d)
      var(y[tmpN[d]:tmpN[d+1], p], na.rm = TRUE)
    ))

  for ( p in 1:P){
    for( d in 1:D){
      temp[d,p] <- (varY_PD[d,p] -sigmaAlpha - sigma2_err[p])/(sigmaBeta)
    }
  }

  out <- rowMeans(temp)
  out[which(out <= 0)] <- 0.2
  return(out)
}

#################################
## full conditionals functions ##
#################################

#-- all nuisance parameters (those of alpha, beta and mix components)--#

variance_fc <- function(param, l_param, nu1, nu2,
                        tau, mu_param, mu_mu_param){
  nu1_star <- (l_param )/2 + nu1
  if( nu1_star <= 1){ nu1_star <- 1}
  if( is.infinite(nu1_star)){ nu1_star <- 100}
  nu2_star <- nu2 + (tau * sum( (param -mu_param)^2, na.rm = T ) +
                       (mu_param - mu_mu_param)^2)/(2 * tau)
  if( nu2_star <= 1){ nu2_star <- 1}
  if( is.infinite(nu2_star)){ nu2_star <- 100}
  out <- 1/rgamma( 1, shape = nu1_star, rate = nu2_star)
  return(out)
}


variance_fc_d <- function(param, l_param, nu1, nu2,
                          tau, mu_param, scale_Fact_d){
  nu1_star <- (l_param )/2 + nu1
  if( is.na(nu1_star)){ nu1_star <- 3}
  if( nu1_star <= 1){ nu1_star <- 1}
  if( is.infinite(nu1_star)){ nu1_star <- 100}

  nu2_star <- nu2 + sum((param - mu_param)^2, na.rm = TRUE)/(2)
  if( is.na(nu2_star)){ nu2_star <- 3}
  if( nu2_star <= 1){ nu2_star <- 1}
  if( is.infinite(nu2_star)){ nu2_star <- 100}

  out <- 1/rgamma( 1, shape = nu1_star, rate = nu2_star)
  OUT <- list( out = out,  nu1_star =  nu1_star, nu2_star = nu2_star )
  return(OUT)
}
variance_fc_di <- function(param, nu1, nu2,
                           mu_param, p_est_d, n){
  nu1_star <- nu1 + sum( p_est_d * n ,na.rm = T)/2
  nu2_star <- nu2 + sum(p_est_d * (param - mu_param)^2, na.rm = TRUE)/(2)
  out <- 1/rgamma( 1, shape = nu1_star, rate = nu2_star)
  return(out)
}



mean_fc <- function( tau, sigma2, param, l_param,
                     mu_mu_param){
  sigma2_star <- (tau * sigma2)/(tau * l_param + 1)
  mu_star <- (tau * sum(param, na.rm = T) + mu_mu_param)/(tau * sigma2)
  mu_star <- sigma2_star * mu_star
  out <- rtruncnorm(1, 0, Inf, mu_star, sqrt(sigma2_star))
  return(out)
}

#--main parameters (alpha, beta)--#

alpha_fc <- function( sigma2_a, beta_p, n, mu_a,
                      z, sigma2_err_p, y_p){
  sigma2_star <- (sigma2_err_p * sigma2_a)/(n * sigma2_a + sigma2_err_p)
  mu_star <- (sum(y_p -  beta_p * z, na.rm = T))/(sigma2_err_p) + mu_a/sigma2_a
  mu_star <- sigma2_star * mu_star
  out <- rtruncnorm(1, 0, Inf, mu_star, sqrt(sigma2_star))
  OUT <- list( out = out,  mu_star =  mu_star, sigma2_star = sigma2_star )
  return(OUT)
}

beta_fc <- function( sigma2_b, alpha_p, n, mu_b,
                     z, sigma2_err_p, y_p){

  sigma2_star <- (sigma2_err_p * sigma2_b)/( sum(z^2, na.rm = T) *
                                               sigma2_b + sigma2_err_p)
  mu_star <- (sum( (y_p -  alpha_p) * z, na.rm = T))/
    (sigma2_err_p) + mu_b/sigma2_b
  mu_star <- sigma2_star * mu_star
  out <- rtruncnorm(1, 0, Inf, mu_star, sqrt(sigma2_star))
  OUT <- list( out = out, mu_star =  mu_star, sigma2_star = sigma2_star )
  return(OUT)
}

#--latent intakes (z) for a given component (d)--#

z_fc <- function( sigma2_d, mu_d, sigma2_err,
                  beta, y_i, alpha, P, scale_Fact_d){

  sigma2_d <- sigma2_d * scale_Fact_d
  sigma2_star <- sum( ((beta^2) * sigma2_d +
                         sigma2_err / P)/(sigma2_d * sigma2_err), na.rm = T )
  sigma2_star <- (sigma2_star)^(-1)
  mu_star <- sum((beta * (y_i - alpha)) / sigma2_err, na.rm = T) +
    (mu_d / sigma2_d)
  mu_star1 <- mu_star
  mu_star <- sigma2_star * mu_star
  sigma2_star <-  sqrt(sigma2_star)

  return(rtruncnorm( 1, 0, Inf, mu_star, sigma2_star))
}

#--error variances (sigma_p)--#

sigma2_err_fc <- function( n, theta1, theta2, y_p,
                           alpha_p, beta_p, z){
  par1 <- (n / 2) + theta1
  par2 <- 0.5 * sum( (y_p - alpha_p -(beta_p * z))^2, na.rm = T ) +theta2
  out <- 1/rgamma(1, shape = par1 , rate = par2)
  OUT <- list( out = out,  nu1_star = par1, nu2_star = par2)
  return(OUT)
}

#--log post MCMC weigths--#
logPost_MCMC_wcum <- function(prob_nc, D, n, labelNew ){ # non cum probs
  out <- sapply( 1:n, function(x) log(prob_nc)[x,labelNew[x]]  )
  out <- sum(out)
  return(out)
}

cauchit_probs <- function(yNew, theta, D){ # Ynew is for a single obs i
  scaling_y <- yNew * theta[-1,1]  # vector length P
  tmp <- sapply(1:(D-1), function(d)
    (1/pi)*( (atan(theta[1,d] + sum(scaling_y, na.rm = T)) + pi/2) )
  ) # length D-1
  tmp2 <- c( tmp[1], diff(tmp) )
  TMP <- c( tmp2, 1 - sum(tmp2,na.rm = T))
  if (sum(c(TMP) < 1) == T ){ TMP <- -Inf}
  return(TMP) # output is a vector of length D
}


dim_range_prob <- function(prob_c, D){
  opz1 <- which(prob_c >=0.5)
  if( length(opz1) != 0 ){
    out <- opz1}else{
      tmp2 <- cbind(seq(1,D), prob_c)
      tmp2 <- tmp2[order(tmp2[,2], decreasing = T),]
      tmp3 <- which(cumsum(tmp2[,2]) >=0.5)[1]
      out <- tmp2[1:tmp3,1]
    }
  return(out)
}

