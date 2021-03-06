% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jc69.2s.R
\name{jc69.rt}
\alias{jc69.rt}
\alias{jc69.rtlnL}
\alias{jc69.rtlnP}
\title{The JC69 model on a pairwise sequence alignment on t and r}
\usage{
jc69.rtlnL(r, t, x, n)

jc69.rtlnP(r, t, x, n, shape.r, rate.r, shape.t, rate.t)
}
\arguments{
\item{r}{numeric, the evolutionary rate}

\item{t}{numeric, the divergence time}

\item{x}{numeric, the number of differences in the alignment}

\item{n}{numeric, the number of sites in the alignment}

\item{shape.r, rate.r, shape.t, rate.t}{numeric, the shape and rate parameters
for the gamma priors on r and t}
}
\description{
The JC69 model on a pairwise sequence alignment on t and r
}
\details{
\code{jc69.rtlnL} and \code{jc69.rtlnP} compute the log-likelihood
  and log-posterior respectively.
}
\examples{
# This reproduces fig. 1 in dos Reis et al. (2016)
data(humanorang)
par(mfrow=c(1,3))
n <- 100
r <- seq(from=0, to=6, len=n) * 1e-3
t <- seq(from=0, to=50, len=n)
rt <- expand.grid(t=t, r=r)
shape.r <- 2; rate.r <- .5 * 1e3  # diffuse prior with mean 4e-3 s/s/My
shape.t <- 16.2; rate.t <- .72    # informative prior with mean 22.5 Ma
# MLE of molecular distance, note d here is distance from tip to root
d.mle <- jc69.2smle(humanorang$x, humanorang$n) / 2

# Joint gamma prior on r and t:
z1 <- dgamma(rt$t, shape.t, rate.t) * dgamma(rt$r, shape.r, rate.r)
z1 <- matrix(z1, ncol=n)
image(t, r, z1, xlab="Root age (Ma)", ylab="Molecular rate (s/s/My)", main="Prior", col=rev(heat.colors(12)), las=1)
contour(t, r, z1, add=TRUE, drawlabels=FALSE)

# Joint likelihood on r and t:
z2 <- exp(jc69.rtlnL(rt$r, rt$t, x=humanorang$x, n=humanorang$n))
z2 <- matrix(z2, ncol=n)
image(t, r, z2, xlab="Root age (Ma)", ylab="Molecular rate (s/s/My)", main="Likelihood", col=rev(heat.colors(12)), las=1)
contour(t, r, z2, add=TRUE, drawlabels=FALSE)
curve(d.mle/x, col="blue", lwd=2, add=TRUE)

# Joint posterior on r and t:
z3 <- exp(jc69.rtlnP(rt$r, rt$t, x=humanorang$x, n=humanorang$n, shape.r, rate.r, shape.t, rate.t))
z3 <- matrix(z3, ncol=n)
image(t, r, z3, xlab="Root age (Ma)", ylab="Molecular rate (s/s/My)", main="Posterior", col=rev(heat.colors(12)), las=1)
contour(t, r, z3, add=TRUE, drawlabels=FALSE)
curve(d.mle/x, col="blue", lwd=2, add=TRUE)

}
\references{
dos Reis M, and Yang Z. (2013) \emph{The unbearable uncertainty of Bayesian
divergence time estimation}. J. Syst. Evol., 51: 30--43.

dos Reis M, Donoghue PCJ, and Yang Z. (2016) \emph{Bayesian molecular clock
dating of species divergences in the genomics era.} Nat. Rev. Genet., 17:
71--80.
}
