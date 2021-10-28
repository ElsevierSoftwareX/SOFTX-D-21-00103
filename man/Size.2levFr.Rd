\name{Size.2levFr}
\alias{Size.2levFr}
\title{Sample Size Calculator for Two-level Fractional Factorial Design}
\description{
This function computes sample size for two-level fractional factorial design to detect a certain standardized effect size with power at the significance level.
The model for fractional factorial design contains only main effects in resolution III and IV.
}
\usage{
Size.2levFr(nfactor, nfraction, interaction = FALSE, delta_type = 1, 
    delta = c(1, 0, 1), alpha = 0.05, beta = 0.2, maxsize = 1000) 
}
\arguments{
\item{nfactor}{the number of factor.}
\item{nfraction}{the number of fraction. For example, when a model is \eqn{2^(k-p)}, k is the number of factor and p is the number of fraction. It is called a \eqn{1/2^p} fraction of the \eqn{2^k} design.} 
\item{interaction}{specifies whether two-way interaction effects are included in a model with the main effects. When \code{interaction = TRUE}, 
two-way interaction effects are include in a model.}
\item{delta_type}{specifies the type of standardized effect size: 1 for standard deviation type and 2 for range type.}
\item{delta}{vector of effect sizes: \code{delta[1]} for main effects, \code{delta[2]} for two-way interaction effects, and \code{delta[3]} for standard deviation of noise.
When \code{interaction=FALSE}, \code{delta[2]} is 0.}
\item{alpha}{Type I error.}
\item{beta}{Type II error.}
\item{maxsize}{tolerance for sample size.}
}
\value{
\item{model}{a character vector expressing a model. The main effects are expressed by the upper-case letters
of the Roman alphabet, and two-way interaction effects are denoted by * operator for pairs of the main effects.}
\item{n}{optimal sample size.}
\item{Delta}{a vector of minimal detectable standardized effect sizes.}
}
\details{
This function computes sample size in two-level fractional factorial design to detect a certain standardized effect size \code{delta} with power \code{1-beta} at the significance level \code{alpha}.
}
\references{
R. V. Lenth (2006-9). Java Applets for Power and Sample Size[Computer software]. Retrieved March 27, 2018 from \url{https://homepage.divms.uiowa.edu/~rlenth/Power/}.

Y. B. Lim (1998). Study on the Size of Minimal Standardized Detectable Difference in Balanced Design of Experiments.
\emph{Journal of the Korean society for Quality Management}, \bold{26(4)}, 239--249.

M. A. Kastenbaum, D. G. Hoel and K. O. Bowman (1970) Sample size requirements : one-way analysis of variance, \emph{Biometrika}, \bold{57(2)}, 421--430.

D. C. Montgomery (2013) Design and analysis of experiments. John Wiley & Sons.
} 
\seealso{
\code{\link{Size.Full}}, \code{\link{Size.Split}}, \code{\link{Size.Block}}.
}
\examples{
# only main effects
model1 <- Size.2levFr(nfactor=3, nfraction=1, interaction=FALSE,
    delta_type=1, delta=c(1, 0, 1), alpha=0.05, beta=0.2)
model1$model
model1$n
model1$Delta

# including two-way interaction effects
model2 <- Size.2levFr(nfactor=5, nfraction=1, interaction=TRUE,
    delta_type=1, delta=c(1, 1, 1), alpha=0.05, beta=0.2)
}

