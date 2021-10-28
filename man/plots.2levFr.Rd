\name{plots.2levFr}
\alias{plots.2levFr}
\title{Diagnosis Graphs for Sample Size of Two-level Fractional Factorial Design}
\description{
This function produces graphs between the sample size, power and the detectable standardized effect size of two-level fractional factorial design.
}
\usage{
plots.2levFr(nfactor, nfraction, interaction = FALSE, delta_type = 1, 
    delta = c(1, 0, 1), deltao = NULL, alpha = 0.05, beta = 0.2, type = 1, 
    maxsize = 1000) 
}
\arguments{
\item{nfactor}{the number of factor.}
\item{nfraction}{the number of fraction. For example, when a model is \eqn{2^(k-p)}, k is the number of factor and p is the number of fraction. It is called a \eqn{1/2^p} fraction of the \eqn{2^k} design.} 
\item{interaction}{specifies whether two-way interaction effects are included in a model with the main effects. When \code{interaction = TRUE}, 
two-way interaction effects are include in a model.}
\item{delta_type}{specifies the type of standardized effect size: 1 for standard deviation type and 2 for range type.}
\item{delta}{vector of effect sizes: \code{delta[1]} for main effects, \code{delta[2]} for two-way interaction effects, and \code{delta[3]} for standard deviation of noise.
When \code{interaction=FALSE}, \code{delta[2]} is 0.}
\item{deltao}{the minimal detectable standardized effect size for power vs the sample size plot when \code{type = 3}.}
\item{alpha}{Type I error.}
\item{beta}{Type II error.}
\item{type}{graph type: 1 for Power vs Delta plot, 2 for Delta vs Sample size plot, and 3 for Power vs Sample size plot.}
\item{maxsize}{tolerance for sample size.}
}
\value{
plot of Power vs Delta, Delta vs Sample size, or Power vs Sample size according to \code{type}.
}
\details{
This function produces graph between the sample size, power 1-\code{beta} and the detectable standardized effect size \code{delta} of two-level fractional factorial design.
According to \code{type}, it displays plot of Power vs Delta, Delta vs Sample size, or Power vs Sample size.
}
\seealso{
\code{\link{plots.Full}}, \code{\link{plots.Split}}, \code{\link{plots.Block}}.
}
\examples{
# plot of Power vs Delta for two-level fractional factorial design
# without the interaction effects
plots.2levFr(nfactor=3, nfraction=1, interaction=FALSE,
    delta_type=1, delta=c(1, 0, 1), alpha=0.05, beta=0.2, type=1)
# plot of Power vs Sample size for two-level fractional factorial design
# with the interaction effects
plots.2levFr(nfactor=5, nfraction=1, interaction=TRUE,
    delta_type=1, delta=c(1, 1, 1), deltao=1, alpha=0.05, beta=0.2, type=3)
}
