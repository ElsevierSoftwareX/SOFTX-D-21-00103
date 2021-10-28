\name{Size.Split}
\alias{Size.Split}
\title{Sample Size Calculator for Split-Plot Design}
\usage{
Size.Split(whole.factor.lev, split.factor.lev, interaction = FALSE, 
    delta_type = 1, delta = c(1, 0, 1, 1), alpha = 0.05, beta = 0.2, 
    maxsize = 1000)
}
\description{
This function computes sample size for split-plot design to detect a certain standardized effect size with power at the significance level.
}
\arguments{
\item{whole.factor.lev}{vector of the numbers of levels for each whole factor.}
\item{split.factor.lev}{vector of the numbers of levels for each split factor.}
\item{interaction}{specifies whether two-way interaction effects are included in a model with the main effects. When \code{interaction = TRUE}, 
two-way interaction effects are include in a model.}
\item{delta_type}{specifies the type of standardized effect size: 1 for standard deviation type and 2 for range type.}
\item{delta}{vector of effect sizes: \code{delta[1]} for main effects, \code{delta[2]} for two-way interaction effects, and \code{delta[3]} and \code{delta[4]} 
 for standard deviation of whole-plot noise and subplot noise, respectively. When \code{interaction=FALSE}, \code{delta[2]} is 0.}
\item{alpha}{Type I error.}
\item{beta}{Type II error.}
\item{maxsize}{tolerance for sample size.}
}
\value{
\item{model}{a character vector expressing a model. The whole factor effects and the split factor effects are expressed by the lower-case letters and sequential upper-case letters of the Roman alphabet, and two-way interaction effects are denoted by * operator for pairs of the those effects.}
\item{n}{optimal sample size.}
\item{Delta}{a vector of minimal detectable standardized effect sizes.}
}
\details{
This function computes sample size in split-plot design to detect a certain standardized effect size \code{delta} with power \code{1-beta} at the significance level \code{alpha}.
The number of whole-plot factors and split plot factors are up to 2 in the current package version.
The linear model for the split-plot design is
\deqn{y_{ijklm} = \mu + \tau_i + \beta_j + \gamma_k + (\beta\tau)_{ik} + \theta_{ijk} + \delta_l + \lambda_m + (\delta\lambda)_{im} + (\beta\delta)_{jl} + 
       (\beta\lambda)_{jm} + (\gamma\delta)_{kl} + (\delta\lambda)_{lm} + \epsilon_{ijklm}}
where \eqn{\tau_i} is the replicate effect, \eqn{\beta_j, \gamma_k} is the whole-plot main effects, \eqn{\theta_{ijk}} is the whole-plot error, 
\eqn{\delta_l, \lambda_m} is the subplot main effects, and \eqn{\epsilon_{ijklm}} is the subplot error.
}
\references{
R. V. Lenth (2006-9). Java Applets for Power and Sample Size[Computer software]. Retrieved March 27, 2018 from \url{https://homepage.divms.uiowa.edu/~rlenth/Power/}.

Y. B. Lim (1998). Study on the Size of Minimal Standardized Detectable Difference in Balanced Design of Experiments.
\emph{Journal of the Korean society for Quality Management}, \bold{26(4)}, 239--249.

M. A. Kastenbaum, D. G. Hoel and K. O. Bowman (1970) Sample size requirements : one-way analysis of variance, \emph{Biometrika}, \bold{57(2)}, 421--430.

D. C. Montgomery (2013) Design and analysis of experiments. John Wiley & Sons.
} 
\seealso{
\code{\link{Size.Full}}, \code{\link{Size.2levFr}}, \code{\link{Size.Block}}.
}
\examples{ 
# only main effects
splitmodel1 <- Size.Split(whole.factor.lev=c(2, 2), split.factor.lev=c(2, 2), interaction=FALSE,
    delta_type=1, delta=c(1, 0, 1, 1), alpha=0.05, beta=0.2)
splitmodel1$model
splitmodel1$n
splitmodel1$Delta

# including two-way interaction effects
splitmodel2 <- Size.Split(whole.factor.lev=c(2, 2), split.factor.lev=c(2, 2), interaction=TRUE,
    delta_type=1, delta=c(1, 1, 1, 1), alpha=0.05, beta=0.2)
splitmodel2
}
