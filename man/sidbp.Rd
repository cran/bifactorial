\name{sidbp}
\alias{sidbp}
\docType{data}
\title{Data of sitting distolic blood pressure (SiDBP)}
\description{
These data have been simulated with respect to the descriptive
statistics given in the bifactorial hypertension clinical trial reported
by Hung (2000). Various distributional properties have been realized for
normal and skewed cases with equal variances as well as with linearly
increasing variances. The latter means that the coefficient of variation
is held constant over the treatment groups. The group defining variable
is named \code{cb} in the data set. It has got the levels
\code{(0,0)}, \code{(1,0)}, \code{...}, \code{(2,3)} according to the respective dose
combinations. The data vectors with the different distributional
properties are \code{ynrmhom}, \code{ynrmhet}, \code{yloghom} and \code{yloghet}.
}
\format{
  A data frame with 738 observations on the following 5 variables.
  \describe{
    \item{\code{cb}}{a factor with levels \code{(0,0)},
      \code{(0,1)}, \code{...}, \code{(2,3)}}
    \item{\code{ynrmhom}}{A vector of normal and homoscedastic data.}
    \item{\code{ynrmhet}}{A vector of normal and heteroscedastic data.}
    \item{\code{yloghom}}{A vector of lognormal and homoscedastic data.}
    \item{\code{yloghet}}{A vector of lognormal and heteroscedastic data.}
  }
}
\references{
  Hung HMJ (2000): Evaluation of a combination drug with multiple doses
  in unbalanced factorial design clinical trials. Statistics in Medicine
  19, pp. 2079-2087
}
\examples{
data(sidbp)
}
\keyword{datasets}
