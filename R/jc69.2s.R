#' The JC69 model on a pairwise sequence alignment
#'
#' @param  d numeric, the molecular distance
#' @param  x numeric, the number of differences in the alignment
#' @param  n numeric, the number of sites in the alignment
#' @param  shape numeric, the shape of the gamma prior on d
#' @param  rate numeric, the rate of the gamma prior on d
#' @param  C numeric, a scaling constant
#'
#' @details \code{jc692s.mle} calulates the MLE of d. \code{jc69.2slnL} computes
#'   the log-likelihood and \code{jc69.2slnP} computes the log-posterior
#'   assuming a gamma prior on \code{d}.
#'
#' @return \code{jc69.2slnL} and \code{jc69.2slnP} return vectors of
#'   log-likelihood and log-posterior values respectively.
#'
#' @references Yang Z (2014) \emph{Molecular evolution: A Statistical approach},
#'   Oxford University Press.
#'
#' @author Mario dos Reis
#'
#' @examples
#' data(humanorang)
#'
#' par(mfrow=c(1,2))
#' # Plot the likelihood:
#' curve(exp(jc69.2slnL(x, x=humanorang$x, n=humanorang$n)),
#'   n=5e2, from=0, to=0.2, ylab="Likelihood", xlab="Molecular distance (d)")
#' abline(v=jc69.2smle(x=humanorang$x, n=humanorang$n), col="red", lty=2)
#'
#' # Plot the posterior:
#' curve(exp(jc69.2slnP(x, x=humanorang$x, n=humanorang$n, shape=1, rate=5)),
#'   n=5e2, from=0, to=0.2, ylab="Posterior", xlab="Molecular distance (d)")
#'
#' # Does not integrate to one:
#' integrate(function(x) exp(jc69.2slnP(x, x=humanorang$x, n=humanorang$n,
#'   shape=1, rate=5)), lower=0, upper=Inf, abs.tol=0)
#'
#' # Integrates to one:
#' integrate(function(x) exp(jc69.2slnP(x, x=humanorang$x, n=humanorang$n,
#'   shape=1, rate=5, C=1/5.167762e-131)), lower=0, upper=Inf, abs.tol=0)
#'
#' @name jc69.2s
NULL

# This is the MLE of the distance
#' @rdname jc69.2s
#' @export
jc69.2smle <- function(x, n) {
  -3/4 * log(1 - 4/3 * x/n)
}

# This is the JC69 log-likelihood
#' @rdname jc69.2s
#' @export
jc69.2slnL <- function(d, x, n)  {
  x*log(3/4 - 3*exp(-4*d/3)/4) + (n - x) *log(1/4 + 3*exp(-4*d/3)/4)
}

# This is the JC69 log-posterior (with gamma prior on d)
#' @rdname jc69.2s
#' @importFrom stats dgamma
#' @export
jc69.2slnP <- function(d, x, n, shape, rate, C=1) {
  jc69.2slnL(d=d, x=x, n=n) + dgamma(x=d, shape=shape, rate=rate, log=TRUE) + log(C)
}

#' The JC69 model on a pairwise sequence alignment on t and r
#'
#' @param r numeric, the evolutionary rate
#' @param t numeric, the divergence time
#' @param x numeric, the number of differences in the alignment
#' @param n numeric, the number of sites in the alignment
#' @param shape.r,rate.r,shape.t,rate.t numeric, the shape and rate parameters
#'   for the gamma priors on r and t
#'
#' @details \code{jc69.rtlnL} and \code{jc69.rtlnP} compute the log-likelihood
#'   and log-posterior respectively.
#'
#' @references
#' dos Reis M, and Yang Z. (2013) \emph{The unbearable uncertainty of Bayesian
#' divergence time estimation}. J. Syst. Evol., 51: 30--43.
#'
#' dos Reis M, Donoghue PCJ, and Yang Z. (2016) \emph{Bayesian molecular clock
#' dating of species divergences in the genomics era.} Nat. Rev. Genet., 17:
#' 71--80.
#'
#' @examples
#' # This reproduces fig. 1 in dos Reis et al. (2016)
#' data(humanorang)
#' par(mfrow=c(1,3))
#' n <- 100
#' r <- seq(from=0, to=6, len=n) * 1e-3
#' t <- seq(from=0, to=50, len=n)
#' rt <- expand.grid(t=t, r=r)
#' shape.r <- 2; rate.r <- .5 * 1e3  # diffuse prior with mean 4e-3 s/s/My
#' shape.t <- 16.2; rate.t <- .72    # informative prior with mean 22.5 Ma
#' # MLE of molecular distance, note d here is distance from tip to root
#' d.mle <- jc69.2smle(humanorang$x, humanorang$n) / 2
#'
#' # Joint gamma prior on r and t:
#' z1 <- dgamma(rt$t, shape.t, rate.t) * dgamma(rt$r, shape.r, rate.r)
#' z1 <- matrix(z1, ncol=n)
#' image(t, r, z1, xlab="Root age (Ma)", ylab="Molecular rate (s/s/My)", main="Prior", col=rev(heat.colors(12)), las=1)
#' contour(t, r, z1, add=TRUE, drawlabels=FALSE)
#'
#' # Joint likelihood on r and t:
#' z2 <- exp(jc69.rtlnL(rt$r, rt$t, x=humanorang$x, n=humanorang$n))
#' z2 <- matrix(z2, ncol=n)
#' image(t, r, z2, xlab="Root age (Ma)", ylab="Molecular rate (s/s/My)", main="Likelihood", col=rev(heat.colors(12)), las=1)
#' contour(t, r, z2, add=TRUE, drawlabels=FALSE)
#' curve(d.mle/x, col="blue", lwd=2, add=TRUE)
#'
#' # Joint posterior on r and t:
#' z3 <- exp(jc69.rtlnP(rt$r, rt$t, x=humanorang$x, n=humanorang$n, shape.r, rate.r, shape.t, rate.t))
#' z3 <- matrix(z3, ncol=n)
#' image(t, r, z3, xlab="Root age (Ma)", ylab="Molecular rate (s/s/My)", main="Posterior", col=rev(heat.colors(12)), las=1)
#' contour(t, r, z3, add=TRUE, drawlabels=FALSE)
#' curve(d.mle/x, col="blue", lwd=2, add=TRUE)
#'
#' @name jc69.rt
NULL

# This is the bi-dimensional JC69 log-likelihood (on r and t)
#' @rdname jc69.rt
#' @export
jc69.rtlnL <- function(r, t, x, n) {
  jc69.2slnL(d=r*t*2, x=x, n=n)
}

# This is the bi-dimensional JC69 log-posterior (on r and t)
#' @rdname jc69.rt
#' @export
jc69.rtlnP <- function(r, t, x, n, shape.r, rate.r, shape.t, rate.t) {
  jc69.rtlnL(r=r, t=t, x=x, n=n) +
    dgamma(r, shape=shape.r, rate=rate.r, log=TRUE) +
    dgamma(t, shape=shape.t, rate=rate.t, log=TRUE)
}
