\name{avemax}
\alias{avetest}
\alias{avetest-class}
\alias{avetest,carpet-method}
\alias{avetest,cube-method}
\alias{avetest,ANY-method}
\alias{maxtest}
\alias{maxtest-class}
\alias{maxtest,carpet-method}
\alias{maxtest,cube-method}
\alias{maxtest,ANY-method}
\alias{show,maxtest-method}
\alias{show,avetest-method}
\alias{summary,maxtest-method}
\alias{summary,avetest-method}
\title{AVE- and MAX-test}
\description{
  Compute global tests for factorial dose-response designs following
  Hung (2000) or by a bootstrap algorithm.
}
\usage{
avetest(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,...)
maxtest(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,...)
}
\arguments{
  \item{C}{An object of class \code{carpet} or \code{cube}.}
  \item{test}{Either \code{"ttest"} or \code{"ztest"} - the test
    statistic for the inferences to be based on. Use \code{"ztest"} for
    binary data applications.}
  \item{method}{The calculation method - use \code{"bootstrap"} for a
    resampling-based approach and \code{"hung"} for calculations using the multivariate normal distribution.}
  \item{nboot}{The number of bootstrap iterations to use.}
  \item{simerror}{Prespecified simulation standard error.}
  \item{\dots}{Any further arguments.}
}
\details{
  When handling with data from factorial clinical trial designs, one
  is often interested in the question whether dose combinations in the trial
  have got a better effect than all of their component drugs, because
  regulatoric requirements demand a contribution to the efficacy by all components. The decision if any of the tested combination drugs has
  got this property can be based on the AVE- or MAX-statistics
  proposed by Hung, Chi and Lipicky (1993). The hypothesis that this is true for
  none of the combinations is rejected if the largest or the average of
  the min-statistics is sufficiently high. The functions
  \code{avetest} and \code{maxtest} calculate the corresponding p-values
  on \code{carpet} or \code{cube} objects with a new bootstrap
  algorithm, which is default, or by the
  multivariate method for unbalanced designs from Hung (2000). A
  resampling-based method is available also for binary data
  applications. The desired simulation accuracy always needs to be specified by the
  number \code{nboot} of simulations to perform or an upper bound
  \code{simerror} for the simulation standard error. If both are
  given, the two constraints will be held simultaneously. Depending on the type of data, the calculations
  can be based on Student's t-test for metric data or the
  Z-statistic for binary applications.
}
\value{
  An object of class \code{avetest} or \code{maxtext}, respectively, with the
  following slots. The slot \code{name} is available for the MAX-test only.
  \item{p}{p-value for the AVE- or MAX-test.}
  \item{stat}{Observed AVE- or MAX-statistic.}
  \item{test}{Type of test statistic which the AVE- or MAX-test was based on.}
  \item{method}{Algorithm used for the calculation.}
  \item{nboot}{Total number of resampling iterations.}
  \item{simerror}{Simulation standard error.}
  \item{name}{Combination group where the maximum of the min-statistics
    was observed.}
  \item{duration}{Total computing duration in seconds.}
  \item{call}{The function call.}
}
\references{
  Frommolt P, Hellmich M (2009): Resampling in multiple-dose factorial
  designs. Biometrical J 51(6), pp. 915-31

  Hellmich M, Lehmacher W (2005): Closure procedures for monotone
  bi-factorial dose-response designs. Biometrics 61, pp. 269-276

  Hung HMJ, Chi GYH, Lipicky RJ (1993): Testing for the existence of a
  desirable dose combination. Biometrics 49, pp. 85-94

  Hung HMJ, Wang SJ (1997): Large-sample tests for binary outcomes in
  fixed-dose combination drug studies. Biometrics 53, pp. 498-503
  
  Hung HMJ (2000): Evaluation of a combination drug with multiple doses
  in unbalanced factorial design clinical trials. Statist Med
  19, pp. 2079-2087
}
\author{
  Peter Frommolt, University of Cologne \email{peter.frommolt@uni-koeln.de}\cr
  \url{http://portal.ccg.uni-koeln.de}
}
\note{
  The performance of the bootstrap-based approach and the method from Hung
  (2000) has been compared and discussed. All algorithms perform very
  conservative if the means in the marginal treatment groups are close
  for the combinations.
  
}
\seealso{
  \code{\link{bifactorial}}, \code{\link{carpet}}, \code{\link{cube}}, \code{\link{mintest}}, \code{\link{margint}}
}
\examples{
#Hypertension example from Hung (2000)
n<-c(75,75,74,48,74,75,74,49,48,50,48,48)
m<-c(0,1.4,2.7,4.6,1.8,2.8,5.7,8.2,2.8,4.5,7.2,10.9)
s<-rep(7.07,12)
x<-list(12)
for(i in 1:12){
  x[[i]]<-rnorm(n[i],mean=0,sd=1)
  x[[i]]<-((x[[i]]-mean(x[[i]]))*(s[i]/sd(x[[i]])))+m[i]
}
hung<-carpet(x,D=c(2,3))
avetest(hung,test="ttest",nboot=20000)
maxtest(hung,test="ttest",nboot=20000)
}
\keyword{htest}
\keyword{multivariate}
\keyword{design}
\keyword{iteration}
