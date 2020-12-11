\name{multiMarker}
\alias{multiMarker}
\alias{print.multiMarker}
\title{A latent variable model to infer food intake from multiple biomarkers.}

\description{
Implements the multiMarker model via an MCMC algorithm.
}

\usage{
multiMarker(y, quantities,
            niter = 10000, burnIn = 3000,
            posteriors = FALSE, sigmaAlpha = 1,
            nuZ1 = NULL, nuZ2 = NULL,
            nuSigmaP1 = NULL, nuSigmaP2 = NULL, sigmaWprior = 0.000001,
            nuBeta1 = 2, nuBeta2 = 3, tauBeta = 0.1)
}


\arguments{
\item{y}{ A matrix of dimension \eqn{(n\times P)}{(n x P)} storing \eqn{P}{P} biomarker measurements on a set of \eqn{n}{n} observations. Missing values (\code{NA}) are allowed.}

\item{quantities}{ A vector of length \eqn{n}{n} storing the food quantities allocated to each of the \eqn{n}{n} observations in the intervention study data. Missing values (\code{NA}) are not allowed.}

\item{niter}{ The number of MCMC iterations. The default value is \code{niter = 10000}.}

\item{burnIn}{ A numerical value, the number of iterations of the chain to be discarded when computing the posterior estimates. The default value is \code{burnIn = 3000}.}

\item{posteriors}{ A logical value indicating if the full parameter chains should also be returned in output. The default value is \code{posteriors = FALSE}.}

\item{sigmaAlpha}{ Intercepts' hyperparameter (\eqn{\sigma_{\alpha^2}}{\sigma_{\alpha^2}}), see details. The default value is \code{sigmaAlpha = 1}.}

\item{nuZ1, nuZ2}{ Are two vectors of length \eqn{D}{D} storing  hyperparameters for the components' variance parameters. The default values are \code{nuZ1 = nuZ2 = NULL}, corresponding to \code{nuZ1 = (D, D-1, ...,1) } and \code{nuZ2 = (D, D,...,D) }.}

\item{nuSigmaP1, nuSigmaP2}{ Scalar hyperparameters for the error's variance parameters. The default values are \code{nuSigmaP1 = nuSigmaP2 = NULL}, corresponding to \code{nuSigmaP1 = 1} and \code{nuSigmaP2 = n}.}

\item{sigmaWprior}{ A scalar corresponding to the components' weights hyperparameter. The default value is \code{sigmaWprior = 0.1} }

\item{nuBeta1, nuBeta2}{ Scalar hyperparameters for the scaling coefficient's variance parameters. The default values are \code{nuBeta1 = 2} and \code{nuBeta2 = 3}.}

\item{tauBeta}{ A scalar factor for the scaling coefficient's variance parameters. The default value is \code{tauBeta = 0.1}.}
}


\value{
An object of class \code{'multiMarker'} containing the following components:
\item{estimates}{ A list with 9 components, storing posterior estimates of medians, standard deviations and \eqn{95\%}{95\%} credible interval lower and upper bounds for the model parameters:
\itemize{
\item{\code{ALPHA_E}}{ is a matrix of dimension \eqn{(4\times P)}{4 x P} storing the posterior estimates of medians (1st row), standard deviations (2nd row) and \eqn{95\%}{95\%} credible interval lower (3rd row) and upper bounds (4th row) for the \eqn{P}{P} intercept parameters, \eqn{(\alpha_1, \dots, \alpha_P)}{(\alpha_1, ..., \alpha_P)}.}
\item{\code{BETA_E}}{ is a matrix of dimension \eqn{(4\times P)}{4 x P} storing the posterior estimates of medians (1st row), standard deviations (2nd row) and \eqn{95\%}{95\%} credible interval lower (3rd row) and upper bounds (4th row) for the \eqn{P}{P} scaling coefficient parameters, \eqn{(\beta_1, \dots, \beta_P)}{(\beta_1, ..., \beta_P)}.}
\item{\code{SigmaErr_E}}{ is a matrix of dimension \eqn{(4\times P)}{4 x P} storing the posterior estimates of medians (1st row), standard deviations (2nd row) and \eqn{95\%}{95\%} credible interval lower (3rd row) and upper bounds (4th row) for the \eqn{P}{P} error variance parameters, \eqn{(\sigma_1^2, \dots, \sigma_P^2)}{(\sigma_1^2, ..., \sigma_P^2)}.}
\item{\code{SigmaD_E}}{ is a matrix of dimension \eqn{(4\times D)}{4 x D} storing the posterior estimates of medians (1st row), standard deviations (2nd row) and \eqn{95\%}{95\%} credible interval lower (3rd row) and upper bounds (4th row) for the \eqn{D}{D} components' variance parameter, \eqn{(\sigma_1^2, \dots, \sigma_D^2)}{(\sigma_1^2, ..., \sigma_D^2)}.}
\item{\code{Z_E}}{ is a matrix of dimension \eqn{(4\times n)}{4 x n} storing the posterior estimates of medians (1st row), standard deviations (2nd row) and \eqn{95\%}{95\%} credible interval lower (3rd row) and upper bounds (4th row) for the \eqn{n}{n} latent intakes, \eqn{(z_1, \dots, z_n)}{(z_1, ..., z_n)}.}
\item{\code{THETA_Est}}{ is an array of \eqn{((P+1)\times (D-1) \times 4)}{((P+1) x (D-1) x 4)} dimensions composed of \eqn{4}{4} \eqn{((P+1)\times (D-1))}{((P+1) x (D-1))} matrices, storing the posterior estimates of medians (1st matrix), standard deviations (2nd matrix) and \eqn{95\%}{95\%} credible interval lower (3rd matrix) and upper bounds (4th matrix) for the components' weights parameters. In each matrix, the first row reports the values for the components' weights intercept parameter, while the other \eqn{P}{P} rows store those of the weights scaling coefficient parameters, \eqn{(\gamma, \theta_1, \dots,\theta_{D-1})}{(\gamma, \theta_1, ...,\theta_{D-1})}. }
\item{\code{sigmaBeta_E}}{ is a vector containing the posterior estimates of medians, standard deviations and \eqn{95\%}{95\%} credible interval lower and upper bounds for the scaling coefficients' variance parameter (\eqn{\sigma_{\beta}^2}{\sigma_{\beta}^2}).}
\item{\code{muAlpha_E}}{ is a vector containing the posterior estimates of medians, standard deviations and \eqn{95\%}{95\%} credible interval lower and upper bounds for the intercepts' mean parameter (\eqn{\mu_{\alpha}}{\mu_{\alpha}}).}
\item{\code{muBeta_E}}{ is a vector containing the posterior estimates of medians, standard deviations and \eqn{95\%}{95\%} credible interval lower and upper bounds for the scaling coefficients' mean parameter (\eqn{\mu_{\beta}}{\mu_{\beta}}).}
\item{\code{varPHp}}{ Estimated error variance parameter values, \eqn{(\nu_{P1}^{*}, \nu_{P2}^{*})}{(\nu_{P1}*, \nu_{P2}*)}, see References.}
}
}
\item{constants}{ A list with 11 components, storing constant model quantities:
\itemize{
\item{\code{nuZ1, nuZ2}}{ are two vectors of length \eqn{D}{D} storing  hyperparameters for the components' variance parameters, see References.}
\item{\code{sigmaAlpha}}{ is a scalar and it corresponds to the variance of the intercept parameters (\eqn{\sigma_{\alpha^2}}{\sigma_{\alpha^2}}).}
\item{\code{nuSigmaP1, nuSigmaP2}}{ are scalar hyperparameters for the error's variance parameters, see References.}
\item{\code{nuBeta1, nuBeta2}}{ are scalar hyperparameters for the scaling coefficient's variance parameters, see References.}
\item{\code{tauBeta}}{ is a scalar factor for the scaling coefficient's variance parameters, see References.}
\item{\code{x_D}}{ is a vector storing the values for the \eqn{D}{D} food quantities. }
\item{\code{P}}{ is a scalar indicating the number of biomarkers in the data.}
\item{\code{D}}{ is a scalar indicating the number of food quantities in the data.}
\item{\code{n}}{ is a scalar indicating the number of observations the data.}
\item{\code{sigmaWprior}}{ is a scalar and it corresponds to the components' weights hyperparameter, see References.}
\item{\code{y_Median}}{ is a vector storing the observed \eqn{P}{P} biomarker median values.}
\item{\code{y_Var}}{ is a vector storing the observed \eqn{P}{P} biomarker variance values.}
}
}
\item{chains}{ If \code{posteriors = TRUE}, a list with posterior distributions of model parameters is returned:
\itemize{
\item{\code{ALPHA_c}}{ is a matrix of dimension \eqn{(niter-burnIn)\times P}{(niter-burnIn) x P} containing the estimated posterior distributions for the intercept parameters, \eqn{(\alpha_1, \dots, \alpha_P)}{(\alpha_1, ..., \alpha_P)}.}
\item{\code{BETA_c}}{ is a matrix of dimension \eqn{(niter-burnIn)\times P}{(niter-burnIn) x P} containing the estimated posterior distributions for the scaling coefficient parameters, \eqn{(\beta_1, \dots, \beta_P)}{(\beta_1, ..., \beta_P)}.}
\item{\code{SigmaErr_c}}{ is a matrix of dimension \eqn{(niter-burnIn)\times P}{(niter-burnIn) x P} containing the estimated posterior distributions for the error variance parameters, \eqn{(\sigma_1^2, \dots, \sigma_P^2)}{(\sigma_1^2, ..., \sigma_P^2)}.}
\item{\code{SigmaD_c}}{ is a matrix of dimension \eqn{(niter-burnIn)\times D}{(niter-burnIn) x D} containing are the estimated posterior distributions for the components' variance parameters, \eqn{(\sigma_1^2, \dots, \sigma_D^2)}{(\sigma_1^2, ..., \sigma_D^2)}.}
\item{\code{Z_c}}{ is a matrix of dimension \eqn{(niter-burnIn)\times n}{(niter-burnIn) x n} containing the estimated posterior distributions for the latent intakes, \eqn{(z_1, \dots, z_n)}{(z_1, ..., z_n)}.}
\item{\code{THETA_c}}{ is an array of \eqn{ (P+1)\times(D-1)\times(niter-burnIn)}{(P+1) x (D-1) x (niter-burnIn)} dimensions containing the estimated posterior distributions for the components' weights parameters. The first one corresponds to that of the weights intercept parameter, while the other \eqn{P}{P} posterior distributions are those of the weights scaling coefficient parameters. In each one of the \eqn{(niter-burnIn)}{(niter-burnIn)} matrices, the first row reports the values for the components' weights intercept parameter, while the other \eqn{P}{P} rows store those of the weights scaling coefficient parameters, \eqn{(\gamma, \theta_1, \dots,\theta_{D-1})}{(\gamma, \theta_1, ...,\theta_{D-1})}.}
\item{\code{sigmaBeta_c}}{ is a vector containing the estimated posterior distribution for the scaling coefficients' variance parameter (\eqn{\sigma_{\beta}^2}{\sigma_{\beta}^2}).}
\item{\code{muAlpha_c}}{ is a vector containing the estimated posterior distribution for the intercepts' mean parameter (\eqn{\mu_{\alpha}}{\mu_{\alpha}}).}
\item{\code{muBeta_c}}{ is a vector containing the estimated posterior distribution for the scaling coefficients' mean parameter (\eqn{\mu_{\beta}}{\mu_{\beta}}).}
\item{\code{weights_info}}{ is a list containing the acceptance probability values for the weights' parameters, \eqn{(\gamma, \theta_1, \dots,\theta_{D-1})}{(\gamma, \theta_1, ...,\theta_{D-1})}.}
}
}
}

\details{ %%%fix
The function facilitates inference of food intake from multiple biomarkers via MCMC, according to the multiMarker model (D'Angelo et al., 2020). The multiMarker model first learns the relationship between the multiple biomarkers and food quantity data from an intervention study and subsequently allows inference on the latent intake when only biomarker data are available.

Consider a biomarker matrix \eqn{\mathbf{Y}}{Y} of dimension \eqn{(n\times P)}{(n x P)}, storing \eqn{P}{P} different biomarker measurements on \eqn{n}{n} independent observations.
The number of food quantities considered in the intervention study is denoted by \eqn{D}{D}, with the corresponding set being \eqn{\mathbf{X}=(X_1, \dots, X_d, \dots, X_D)}{X=(X_1, ..., X_d, ..., X_D)} and \eqn{X_d < X_{d+1}}{X_d < X_{d+1}}.

We assume that the biomarker measurements are related to an unobserved, continuous intake value, leading to the following factor analytic model:

\deqn{
y_{ip} = \alpha_p + \beta_p z_i +\epsilon_{ip}, \quad \forall \quad  i=1,\dots,n, \quad p = 1, \dots,P,}{
y_{ip} = \alpha_p + \beta_p z_i +\epsilon_{ip}, for  i=1,...,n, and p = 1, ...,P,
}

where the latent variable \eqn{z_i}{z_i} denotes the latent intake of observation \eqn{i}{i}, with \eqn{\mathbf{z}=(z_1, \dots, z_i, \dots, z_n)}{(z_1, ..., z_i, .., z_n)}.
The \eqn{\alpha_p}{\alpha_p} and \eqn{\beta_p}{\beta_p} parameters characterize, respectively, the intercept and the scaling effect for biomarker \eqn{p}{p}. We assume that these parameters are distributed a priori according to 0-truncated Gaussian distributions, with parameters \eqn{(\mu_{\alpha}, \sigma_{\alpha}^2)}{(\mu_{\alpha}, \sigma_{\alpha}^2 )} and \eqn{(\mu_{\beta}, \sigma_{\beta}^2 )}{(\mu_{\beta}, \sigma_{\beta}^2 )} respectively. The error term \eqn{\epsilon_p}{\epsilon_p} is the variability associated with biomarker \eqn{p}{p}. We assume that these errors are normally distributed with 0 mean and variance \eqn{\sigma_p^2}{\sigma_p^2}, which serves as a proxy for the precision of the biomarker.

A mixture of \eqn{D}{D} 0-truncated Gaussian distributions is assumed as prior distribution for the latent intakes. Components are centered around food quantity values \eqn{X_d}{X_d}, and component-specific variances \eqn{\theta_d^2}{\theta_d^2} model food quantity-specific intake variability, with lower values suggesting higher consumption-compliance. Mixture weights are observation-specific and denoted with \eqn{\pi_i = (\pi_{i1}, \dots,\pi_{iD} ) }{\pi_i =(\pi_{i1}, ..., \pi_{iD})}. Given the inherent ordering of the food quantities in the intervention study, an ordinal regression model with Cauchit link function is employed to model the observation-specific weights.

A Bayesian hierarchical framework is employed for the modelling process, allowing quantification of the uncertainty in intake estimation, and flexibility in adapting to different biomarker data distributions. The framework is implemented through a Metropolis within Gibbs Markov chain Monte Carlo (MCMC) algorithm. Hyperprior distributions are assumed on the prior parameters with the corresponding hyperparameter values fixed based on the data at hand, following an empirical Bayes approach.

For more details on the estimation of the multiMarker model, see D'Angelo et al. (2020).
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

#-- Visualize the data --#
par(mfrow= c(2,2))
boxplot(Y[,1] ~ quantities, xlab = "Food quantity", ylab = "Biomarker 1")
boxplot(Y[,2] ~ quantities, xlab = "Food quantity", ylab = "Biomarker 2")
boxplot(Y[,3] ~ quantities, xlab = "Food quantity", ylab = "Biomarker 3")

#-- Fit the multiMarker model --#
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

par(oldpar)
}

\references{
D'Angelo, S. and Brennan, L. and Gormley, I.C. (2020). Inferring food intake from multiple biomarkers using a latent variable model. \emph{arxiv}
}

