\name{plots.Block}
\alias{plots.Block}
\title{Diagnosis Graphs for the number of Blocks of Randomized Complete Block Design}  
\description{
This function produces graphs between the sample size, power and the detectable standardized effect size of randomized complete block design.
}
\usage{
plots.Block(factor.lev, interaction = FALSE, delta_type = 1, delta = c(1, 0, 1), 
    deltao = NULL, alpha = 0.05, beta = 0.2, type = 1, maxsize = 1000) 
}
\arguments{
\item{factor.lev}{vector of the numbers of levels for each factor.}
\item{interaction}{specifies whether two-way interaction effects are included in a model with the main effects. When \code{interaction = TRUE}, 
two-way interaction effects are include in a model.}
\item{delta_type}{specifies the type of standardized effect size: 1 for standard deviation type and 2 for range type.}
\item{delta}{vector of effect sizes: \code{delta[1]} for main effects, \code{delta[2]} for two-way interaction effects, and \code{delta[3]} for standard deviation of noise.
When \code{interaction=FALSE}, \code{delta[2]} is 0.}
\item{deltao}{the minimal detectable standardized effect size for power vs the number of blocks plot when \code{type = 3}.}
\item{alpha}{Type I error.}
\item{beta}{Type II error.}
\item{type}{graph type: 1 for Power vs Delta plot, 2 for Delta vs the Number of Blocks plot, and 3 for Power vs the Number of Blocks plot.}
\item{maxsize}{tolerance for the number of blocks.}
}
\value{
plot of Power vs Delta, Delta vs Number of Blocks, or Power vs Number of Blocks according to \code{type}.
}
\details{
In a randomized complete block design (without replications), the optimal number of blocks need to be determined.
This function produces graph between Number of Block, power 1-\code{beta} and the detectable standardized effect size \code{delta} of randomized complete block design.
According to \code{type}, it displays plot of Power vs Delta, Delta vs Number of Blocks, or Power vs Number of Blocks.
}
\seealso{
\code{\link{plots.Full}}, \code{\link{plots.2levFr}}, \code{\link{plots.Split}}.
}
\examples{
# plot of Power vs Delta for randomized complete block design 
# with 2 factors without the interaction effects
plots.Block(factor.lev=c(2, 2), interaction=FALSE,
    delta_type=1, delta=c(1, 0, 1), alpha=0.05, beta=0.2, type=1)
    
# plot of Power vs Number of Blocks for randomized complete block design 
# with 2 factors with the interaction effects
plots.Block(factor.lev=c(2, 3), interaction=TRUE,
    delta_type=1, delta=c(1, 1, 1), deltao=1.5, alpha=0.05, beta=0.2, type=3)
}