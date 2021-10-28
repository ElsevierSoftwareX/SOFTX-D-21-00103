\name{plots.Split}
\alias{plots.Split}
\title{Diagnosis Graphs for Sample Size of Split-Plot Design}
\description{
This function produces graphs between the sample size, power and the detectable standardized effect size of split-plot design.}
\usage{
plots.Split(whole.factor.lev, split.factor.lev, interaction = FALSE, 
    delta_type = 1, delta = c(1, 0, 1, 1), deltao = NULL, alpha = 0.05, beta = 0.2, 
    type = 1, maxsize = 1000)
}
\arguments{
\item{whole.factor.lev}{vector of the numbers of levels for each whole factor.}
\item{split.factor.lev}{vector of the numbers of levels for each split factor.}
\item{interaction}{specifies whether two-way interaction effects are included in a model with the main effects. When \code{interaction = TRUE}, 
two-way interaction effects are include in a model.}
\item{delta_type}{specifies the type of standardized effect size: 1 for standard deviation type and 2 for range type.}
\item{delta}{vector of effect sizes: \code{delta[1]} for main effects, \code{delta[2]} for two-way interaction effects, and \code{delta[3]} and \code{delta[4]} 
 for standard deviation of whole-plot noise and subplot noise, respectively. When \code{interaction=FALSE}, \code{delta[2]} is 0.}
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
This function produces graph between the sample size, power 1-\code{beta} and the detectable standardized effect size \code{delta} of split-plot design.
According to \code{type}, it displays plot of Power vs Delta, Delta vs Sample size, or Power vs Sample size.
The number of whole-plot factors and split plot factors are up to 2 in the current package version.
}
\seealso{
\code{\link{plots.Full}}, \code{\link{plots.2levFr}}, \code{\link{plots.Block}}.
}
\examples{
# plot of Power vs Delta for split-plot design 
# without the interaction effects
plots.Split(whole.factor.lev=2, split.factor.lev=2, interaction=FALSE,
    delta_type=1, delta=c(1, 0, 1, 1), alpha=0.05, beta=0.2, type=1)
  
# plot of Power vs Sample size for split-plot design 
# with the interaction effects
plots.Split(whole.factor.lev=2, split.factor.lev=2, interaction=TRUE,
    delta_type=1, delta=c(1, 1, 1, 1), deltao=1, alpha=0.05, beta=0.2, type=3)
}