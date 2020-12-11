\name{predict.multiMarker}
\alias{predict.multiMarker}
\title{A latent variable model to infer food intake from multiple biomarker data alone.}

\description{
Implements the multiMarker model via an MCMC algorithm.
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
\item{inferred_E}{ a list with 2 components, storing estimates of medians, standard deviations and \eqn{95\%}{95\%} credible interval lower and upper bounds for:
\itemize{
\item{\code{inferred_intakes}}{ is a matrix of dimension \eqn{(4\times n^*)}{(4 x n*)}, storing the estimates of medians (1st row), standard deviations (2nd row) and \eqn{95\%}{95\%} credible interval lower (3rd row) and upper bounds (4th row) from the conditional distribution of the \eqn{n^{*}}{n*} latent intakes, \eqn{(z_1^{*}, \dots, z_{n^*})}{(z_1*, ..., z_n*)}.}
\item{\code{inferred_Prob}}{ is an array of dimension \eqn{(n^{*}\times D\times 4)}{(n* x D x 4)}, storing estimated median (1st matrix), standard deviation (2nd matrix) and \eqn{95\%}{95\%} credible interval lower (3rd matrix) and upper bound (4th matrix) values for the food quantity probabilities, for each one of the new \eqn{ {n^*}}{n*} observations.}
}
}
\item{chains}{ If \code{posteriors = TRUE}, it contains a list with conditional distributions for:
\itemize{
\item{\code{ZINF}}{ is a matrix of dimension \eqn{n^{*}\times niter}{n* x niter} containing samples from the conditional distributions of the latent intakes, \eqn{(z_1^{*}, \dots, z_{n^*})}{(z_1*, ..., z_n*)}.}
\item{\code{PROBS}}{ is an array of \eqn{ n^{*}\times D \times niter }{n^{*} x D x niter} dimensions containing samples from the conditional distribution for food quantity probabilities, for each observation and food quantity.}
}
}
}

\details{ %%%fix
The function facilitates inference on food intake from multiple biomarkers alone via MCMC, according to the multiMarker model (D'Angelo et al., 2020).

A Bayesian framework is employed for the modelling process, allowing quantification of the uncertainty associated with inferred intake. The framework is implemented through an MCMC algorithm.

For more details, see D'Angelo et al. (2020).
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

#-- Infer intakes from biomarker only data --#
# Number of iterations (and burnIn) set small for example.
infM <- predict(modM, y = YNew, niter = 100, burnIn = 30,
                 posteriors = TRUE)
# niter and burnIn values are low only for example purpose

#-- Extract summary statistics for a given intake --#
obs_j <- 2 # choose which observation to look at
infM$inferred_E$inferred_intakes[, obs_j] #inferred median, standard deviation,
# 0.025 and 0.975 quantiles for the intake of observation obs_j

#-- Example of plot --#
par(mfrow = c(1,1))
hist(infM$chains$ZINF[obs_j, ], breaks = 50,
    ylab = "Density", xlab = "Intake",
    main = "Intake's conditional distribution",
    cex.main = 0.7,
    freq = FALSE) # Inferred condtional distribution of intake for observation obs_j
abline( v = infM$inferred_E$inferred_intakes[1,obs_j], col = "darkred",
lwd = 2 ) # median value
abline( v = infM$inferred_E$inferred_intakes[3,obs_j], col = "grey",
lwd = 2 )
abline( v = infM$inferred_E$inferred_intakes[4,obs_j], col = "grey",
lwd = 2 )
legend( x = "topleft", fill = c("grey", "darkred"), title = "quantiles:",
legend = c("(0.025, 0.975)", "0.5"), bty = "n", cex = 0.7)

mtext(paste("Observation", obs_j, sep = " "), outer = TRUE, cex = 1.5)
par(oldpar)

}

\references{
D'Angelo, S. and Brennan, L. and Gormley, I.C. (2020). Inferring food intake from multiple biomarkers using a latent variable model. \href{https://arxiv.org/abs/2006.02995}{arXiv}.
}

