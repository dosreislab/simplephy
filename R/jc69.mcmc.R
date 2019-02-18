#' MCMC sampling of r and t under JC69
#'
#' @examples
#' data(humanorang)
#' n <- 100
#' r <- seq(from=0, to=6, len=n) * 1e-3
#' t <- seq(from=0, to=50, len=n)
#' rt <- expand.grid(t=t, r=r)
#' shape.r <- 2; rate.r <- .5 * 1e3 # diffuse prior with mean 4e-3 substitutions per My
#' shape.t <- 16.2; rate.t <- .72   # informative prior age of root with mean 22.5 Ma
#' # MLE of molecular distance, note d here is distance from tip to root
#' d.mle <- jc69.2smle(humanorang$x, humanorang$n) / 2
#'
#' # Joint posterior on r and t:
#' z3 <- exp(jc69.rtlnP(rt$r, rt$t, x=humanorang$x, n=humanorang$n, shape.r, rate.r, shape.t, rate.t))
#' z3 <- matrix(z3, ncol=n)
#' image(t, r, z3, xlab="Root age (Ma)", ylab="Molecular rate (s/s/My)", main="Posterior", col=rev(heat.colors(12)), las=1)
#' contour(t, r, z3, add=TRUE, drawlabels=FALSE)
#' curve(d.mle/x/1e-3, col="blue", lwd=2, add=TRUE)
#'
#' # MCMC sampling
#' rt.mcmc <- jc69.rtmcmc(init.r=1e-3, init.t=10, N=1e4, x=humanorang$x, n=humanorang$n,
#'   shape.r, rate.r, shape.t, rate.t, w.r=2.5e-3, w.t=20)
#' # Add trace to contour plot:
#' points(rt.mcmc$t, rt.mcmc$r, col="blue", pch='.')
#'
#' @export
jc69.rtmcmc <- function(init.r, init.t, N, x, n, shape.r, rate.r, shape.t, rate.t, w.r, w.t) {
  # init.r and init.t are the initial states
  # N is the number of 'generations' the algorithm is run for.
  # w.r and w.t are the step sizes of the proposal densities.
  r <- t <- numeric(N+1)
  r[1] <- init.r; t[1] <- init.t
  rnow <- init.r; tnow <- init.t
  ulnP <- jc69.rtlnP(rnow, tnow, x, n, shape.r, rate.r, shape.t, rate.t)
  ap.r <- ap.t <- 0 # acceptance proportions

  for (i in 1:N) {
    # Propose and accept or reject new r:
    rnew <- abs(rnow + runif(1, -w.r/2, w.r/2))
    ulnPprop <- jc69.rtlnP(rnew, tnow, x, n, shape.r, rate.r, shape.t, rate.t)
    lnalpha <- ulnPprop - ulnP

    if (lnalpha > 0 || runif(1) < exp(lnalpha)) {
      rnow <- rnew; ulnP <- ulnPprop
      ap.r  <- ap.r + 1
    }

    # Propose and accept or reject new t:
    tnew <- abs(tnow + runif(1, -w.t/2, w.t/2))
    ulnPprop <- jc69.rtlnP(rnow, tnew, x, n, shape.r, rate.r, shape.t, rate.t)
    lnalpha <- ulnPprop - ulnP

    if (lnalpha > 0 || runif(1) < exp(lnalpha)) {
      tnow <- tnew; ulnP <- ulnPprop
      ap.t  <- ap.t + 1
    }
    r[i+1] <- rnow;  t[i+1] <- tnow;   # take the sample
  }
  # print out the acceptance proportions
  print(c(ap.r/N, ap.t/N))
  return (list(r=r, t=t, ap.r=ap.r/N, ap.t=ap.t/N)) # return mcmc states
}
