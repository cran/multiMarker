\name{predict.multiMarker}
\alias{predict.multiMarker}
\title{A latent variable framework to predict food intake from multiple biomarkers.}

\description{
Implements the multiMarker framework via an MCMC algorithm.
}

\usage{
\method{predict}{multiMarker}( object, y,
         niter = 10000, burnIn = 3000,
         posteriors = FALSE, ...)
}


\arguments{
\item{object}{ An object of class inheriting from \code{'multiMarker'}.}

\item{y}{ A matrix of dimension \eqn{(n^{*} \times P)}{(n* x P)} storing \eqn{P}{P} biomarker measurements on a set of \eqn{n^{*}}{n*} observations. Missing values (\code{NA}) are allowed.}

\item{niter}{ The number of MCMC iterations. The default value is \code{niter = 10000}.}

\item{burnIn}{ A numerical value, the number of iterations of the chain to be discarded when computing the posterior estimates. The default value is \code{burnIn = 3000}.}

\item{posteriors}{ A logical value indicating if the full parameter chains should also be returned in output. The default value is \code{posteriors = FALSE}.}

\item{...}{ Further arguments passed to or from other methods.}
}


\value{A list with 2 components:
\item{predictions}{ a list with 2 components, storing posterior predictive estimates of medians, standard deviations and \eqn{95\%}{95\%} credible interval lower and upper bounds for:
\itemize{
\item{\code{ZPred_P}}{ is a matrix of dimension \eqn{(4\times n^*)}{(4 x n*)}, storing the posterior predictive estimates of medians (1st row), standard deviations (2nd row) and \eqn{95\%}{95\%} credible interval lower (3rd row) and upper bounds (4th row) for the \eqn{n^{*}}{n*} latent intakes, \eqn{(z_1^{*}, \dots, z_{n^*})}{(z_1*, ..., z_n*)}.}
\item{\code{Prob_P}}{ is an array of dimension \eqn{(n^{*}\times D\times 4)}{(n* x D x 4)}, storing estimated median (1st matrix), standard deviation (2nd matrix) and \eqn{95\%}{95\%} credible interval lower (3rd matrix) and upper bound (4th matrix) values for the posterior predictive food quantity probabilities, for each one of the new \eqn{ {n^*}}{n*} observations.}
}
}
\item{chains}{ If \code{posteriors = TRUE}, it contains a list with posterior predictive distributions for:
\itemize{
\item{\code{ZPRED}}{ is a matrix of dimension \eqn{n^{*}\times niter}{n* x niter} containing samples from the posterior predictive distributions of the latent intakes, \eqn{(z_1^{*}, \dots, z_{n^*})}{(z_1*, ..., z_n*)}.}
\item{\code{PROBS}}{ is an array of \eqn{ n^{*}\times D \times niter }{n^{*} x D x niter} dimensions containing samples from the posterior predictive distribution for food quantity probabilities, for each observation and food quantity.}
}
}
}

\details{ %%%fix
The function facilitates food intake prediction from multiple biomarkers via MCMC, according to the multiMarker framework (D'Angelo et al., 2020).

When new biomarker data are available, the multiMarker framework can use such data to perform intake prediction. Indeed, consider \eqn{n^{*}}{n*} new observations, \eqn{j=1,\dots,n^{*}}{j=1,...,n*}, for which only biomarker data are observed (\eqn{Y^{*}}{Y*}). Here, predicted intake values for new observation \eqn{j}{j}, \eqn{z_j^{*}}{z_j*}, can be sampled from its predictive posterior distribution:

\deqn{p(z_j^{*}\mid \mathbf{z}, y_{j}^{*}) \propto\int p(z_j^{*} \mid \Omega, y_{j}^{*}) p(\Omega \mid \mathbf{z})  \,d\Omega}{
p(z_j* | z, y_j*) ∝ ∫ p(z_j* | \Omega, y_j*) p(\Omega | z)  d\Omega
}

where \eqn{\Omega = \{  \theta, \gamma, c_j, \mu_{z1}, \dots, \mu_{zD}, \sigma_{zq}^2, \dots , \sigma_{zD}^2\}}{\Omega = \{  \theta, \gamma, c_j, \mu_{z1}, \dots, \mu_{zD}, \sigma_{zq}^2, \dots , \sigma_{zD}^2\}} is the set of prediction parameters, and

\deqn{
p(z_j^{*} \mid \Omega, y_{j}^{*}) =\sum_{d=1}^D \pi_d( \theta_d, \gamma, y_j^{*}, c_{jd}) \mathcal{N}_{[0, \infty]} (\mu_{zd} , \sigma_{zd}^2)
}{
p(z_j* | \Omega, y_j*) =\sum_{d=1}^D \pi_d( \theta_d, \gamma, y_j*, c_{jd}) N_[0, \infty] (\mu_{zd} , \sigma_{zd}^2)
}

is the sampling distribution for \eqn{z_j^{*}}{z_j*}, with \eqn{\mu_{zd} = \frac{\mu_z\sigma_d^2\tau_d + X_d\sigma_z^2}{\sigma_d^2\tau_d +\sigma_z^2}}{\mu_{zd} = (\mu_z\sigma_d^2\tau_d + X_d\sigma_z^2)/(\sigma_d^2\tau_d +\sigma_z^2)} and \eqn{\sigma_{zd}^2 =\bigl( \frac{1}{\sigma_d^2\tau_d } + \frac{1}{\sigma_z^2} \bigr)^{-1}}{\sigma_{zd}^2 = ( 1/(\sigma_d^2\tau_d) + 1/\sigma_z^2 )^{-1}}.

A Bayesian framework is employed for the modelling process, allowing quantification of the uncertainty in intake predictions. The framework is implemented through an MCMC algorithm.

For more details on how to perform latent intake predictions, see D'Angelo et al. (2020).
}


\examples{

library(truncnorm)
oldpar <- par(no.readonly =TRUE)

#-- Simulate intervention study biomarker and food quantity data --#

P <- D <- 3; n <- 50
alpha <- rtruncnorm(P, 0, Inf, 4, 1)
beta <- rtruncnorm(P, 0, Inf, 0.001, 0.1)
x <- c(50, 100, 150)
labels_z <- sample(c(1,2,3), n, replace = TRUE)
quantities <- x[labels_z]
sigma_d <- 8
z <- rtruncnorm(n, 0, Inf, x[labels_z], sigma_d)
Y <- sapply( 1:P, function(p) sapply( 1:n, function(i)
  max(0, alpha[p] + beta[p]*z[i] + rnorm( 1, 0, 5) ) ) )

#-- Simulate Biomarker data only --#
nNew <- 20
labels_zNew <- sample(c(1,2,3), nNew, replace = TRUE)
zNew <- rtruncnorm(nNew, 0, Inf, x[labels_zNew], sigma_d)
YNew <- sapply( 1:P, function(p) sapply( 1:nNew, function(i)
  max(0, alpha[p] + beta[p]*zNew[i] + rnorm( 1, 0, 5) ) ) )

#-- Fit the multiMarker model to the intervention study data --#
# Number of iterations (and burnIn) set small for example.
modM <- multiMarker(y = Y, quantities = quantities,
                    niter = 100, burnIn = 30,
                    posteriors = TRUE)
                    # niter and burnIn values are low only for example purposes

#-- Extract summary statistics for model parameters --#
modM$estimates$ALPHA_E[,3] #estimated median, standard deviation,
# 0.025 and 0.975 quantiles for the third intercept parameter (alpha_3)

modM$estimates$BETA_E[,2] #estimated median, standard deviation,
# 0.025 and 0.975 quantiles for the second scaling parameter (beta_2)

#-- Examine behaviour of MCMC chains --#
par(mfrow= c(2,1))
plot(modM$chains$ALPHA_c[,3], type = "l",
xlab = "Iteration (after burnin)", ylab = expression(alpha[3]) )
abline( h = mean(modM$chains$ALPHA_c[,3]), lwd = 2, col = "darkred")

plot(modM$chains$BETA_c[,2], type = "l",
xlab = "Iteration (after burnin)", ylab = expression(beta[2]) )
abline( h = mean(modM$chains$BETA_c[,2]), lwd = 2, col = "darkred")

# compute Effective Sample Size
# library(LaplacesDemon)
# ESS(modM$chains$ALPHA_c[,3]) # effective sample size for alpha_3 MCMC chain
# ESS(modM$chains$BETA_c[,2]) # effective sample size for beta_2 MCMC chain

#-- Predict intakes from biomarker only data --#
# Number of iterations (and burnIn) set small for example.
predM <- predict(modM, y = YNew, niter = 100, burnIn = 30,
                 posteriors = TRUE)
# niter and burnIn values are low only for example purpose

#-- Extract summary statistics for a given intake --#
obs_j <- 2 # choose which observation to look at
predM$predictions$ZPred_P[, obs_j] #predicted median, standard deviation,
# 0.025 and 0.975 quantiles for the intake of observation obs_j

#-- Example of posterior predictive plots --#
obs_j <- 2 # choose which observation to look at
par(mfrow=c(1,2), oma = c(0, 0, 2, 0))
barplot(predM$predictions$Prob_P[obs_j,,1],
        xlab = "Food quantity", ylab = "Probability",
        main = "Median posterior probability \nof consuming intervention study \nfood quantities",
cex.main = 0.7) #  Posterior food quantity probabilities
axis(1, at = seq(1, D), labels = x)

hist(predM$chains$ZPRED[obs_j, ], breaks = 50,
    ylab = "Density", xlab = "Intake",
    main = "Posterior predictive distribution \nof intake",
    cex.main = 0.7,
    freq = FALSE) # Posterior predictive distribution of intake for observation obs_j
abline( v = predM$predictions$ZPred_P[1,obs_j], col = "darkred",
lwd = 2 ) # median value
abline( v = predM$predictions$ZPred_P[3,obs_j], col = "grey",
lwd = 2 )
abline( v = predM$predictions$ZPred_P[4,obs_j], col = "grey",
lwd = 2 )
legend( x = "topleft", fill = c("grey", "darkred"), title = "quantiles:",
legend = c("(0.025, 0.975)", "0.5"), bty = "n", cex = 0.7)

mtext(paste("Observation", obs_j, sep = " "), outer = TRUE, cex = 1.5)
par(oldpar)

}

\references{
D'Angelo, S. and Brennan, L. and Gormley, I.C. (2020). Inferring food intake from multiple biomarkers using a latent variable model. \href{https://arxiv.org/abs/2006.02995}{arXiv}.
}

