\name{multinf}
\alias{mintest}
\alias{mintest-class}
\alias{mintest,carpet-method}
\alias{mintest,cube-method}
\alias{mintest,ANY-method}
\alias{show,mintest-method}
\alias{summary,mintest-method}
\alias{plot,mintest,missing-method}
\alias{margint}
\alias{margint-class}
\alias{margint,carpet-method}
\alias{margint,cube-method}
\alias{margint,ANY-method}
\alias{show,margint-method}
\alias{summary,margint-method}
\alias{plot,margint,missing-method}
\title{Multiple inference}
\description{
  Compute adjusted p-values and simultaneous confidence intervals for
  given bi- and trifactorial design data.
}
\usage{
mintest(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,...)
margint(C,test=NULL,method="bootstrap",nboot=NULL,simerror=NULL,
        alpha=0.05,...)
}
\arguments{
  \item{C}{An object of class \code{carpet} or \code{cube}.}
  \item{method}{The calculation method - use \code{"bootstrap"} for a
    resampling-based approach, \code{"hung"} for the min-test approach
    of Hung (2000) and \code{"tdistr"} for interval calculations based
    on the multivariate t-distribution.}
  \item{test}{Either \code{"ttest"} or \code{"ztest"} - the test
    statistic for the inferences to be based on. Use \code{"ztest"} for
    binary data applications.}
  \item{alpha}{Simultaneous level of the confidence intervals.}
  \item{nboot}{Number of resampling iterations to use.}
  \item{simerror}{Prespecified simulation standard error.}
  \item{\dots}{Any further arguments.}
}
\details{
  The generic functions \code{mintest} and \code{margint} calculate
  adjusted p-values and simultaneous confidence intervals for the test
  of parametric differences between prespecified treatment groups on bi-
  or trifactorial design clinical trials. If an object of class
  \code{carpet} is commited, \code{mintest} will return adjusted
  p-values for the min-test on combination superiority in bifactorial
  clinical trial designs (Laska and Meisner, 1989). The alternative
  hypothesis of this test is that the detected effect size for the
  combination treatment is better than for both single component groups;
  i.e. the test results in only one p-value for each combination. The
  generic function \code{margint} will, when applied to \code{carpet}
  objects, return simultaneous confidence intervals for the parametric
  differences between each combination treatment group and its
  respective components. Depending on the type of data, the calculations
  can be based on Student's t-test for metric data or the Z-statistic
  for binary applications.

  By default, the calculations are performed by a resampling-based
  approach. The desired simulation accuracy always needs to be specified
  by the number \code{nboot} of bootstrap iterations to perform or an
  upper bound \code{simerror} for the simulation standard error. If both
  are given, the two constraints will be held simultaneously. On the
  other hand, the multivariate normal approach for unbalanced designs
  from Hung (2000) is available when the argument \code{method} is set
  to the value \code{"hung"}. For the trifactorial case, no such
  approach is available and thus the calculations are based on the
  bootstrap approach, performing a generalized min-test on the data, if
  an object of class \code{cube} is commited. The interval calculations
  are based on the multivariate t-distribution if \code{"tdistr"} is
  specified.

  In the classical approach to the min-test, a normality assumption for
  the data is used and the desired critical values are calculated using
  quantiles of the multivariate t-distribution. However, this method
  fails when handling with data that are skewed or heteroscedastic over
  the treatment groups. When using the bootstrap, only the empirical
  distribution of the data is used and thus the results are always
  valid, provided that a sufficiently large samples are available. When
  handling with data from bifactorial clinical trial designs, bootstrap
  methods need much less analytical framework on the distributional
  properties of the tests than if the approach given by Hung (2000) is
  used. In particular, the restriction to only two compounds is not
  needed and binary data applications can be handled analogously. The
  theory of resampling-based multiple testing has been extensively
  discussed by Westfall and Young (1993).
  
  The calculation of simultaneous confidence intervals is much easier
  because the c.d.f. of the min-statistic is not needed. Hence this is
  leading to an ordinary multiple contrast problem.
}
\value{
  An object of class \code{mintest} or \code{margint} with the following slots.
  \item{p}{Adjusted p-values for the respective combination groups.}
  \item{stat}{The observed values of the min-statistics.}
  \item{kiu}{The lower limits of the confidence intervals.}
  \item{kio}{The upper limits of the confidence intervals.}
  \item{alpha}{One minus the nominal coverage probability of the confidence intervals.}
  \item{gnames}{Names of the combination groups.}
  \item{cnames}{The names of the contrasts for comparisons of the
    combinations with their respective components.}
  \item{test}{Type of test statistics that the min-tests were based on.}
  \item{method}{The method used for calculation.}
  \item{nboot}{Number of bootstrap replications used.}
  \item{simerror}{Maximum of the simulation standard errors in the combination groups.}
  \item{duration}{Total computing duration in seconds.}
  \item{call}{Function call.}
}
\references{
  Frommolt P, Hellmich M (2009): Resampling in multiple-dose factorial
  designs. Biometrical J 51(6), pp. 915-31

  Hung HMJ, Chi GYH, Lipicky RJ (1993): Testing for the existence of a
  desirable dose combination. Biometrics 49, pp. 85-94

  Hung HMJ, Wang SJ (1997): Large-sample tests for binary outcomes in
  fixed-dose combination drug studies. Biometrics 53, pp. 498-503

  Hung HMJ (2000): Evaluation of a combination drug with multiple doses in
  unbalanced factorial design clinical trials. Statist Med 19, pp. 2079-2087

  Hellmich M, Lehmacher W (2005): Closure procedures for monotone bi-factorial
  dose-response designs. Biometrics 61, pp. 269-276

  Laska EM, Meisner MJ (1989): Testing whether an identified treatment is
  best. Biometrics 45, pp. 1139-1151

  Snapinn SM (1987): Evaluating the efficacy of a combination
  therapy. Statist Med 6, pp. 657-665
  
  Westfall PH, Young SS (1993): Resampling-based multiple testing. John
  Wiley & Sons, Inc., New York
}
\author{Peter Frommolt, University of Cologne \email{peter.frommolt@uni-koeln.de}\cr
  \url{http://portal.ccg.uni-koeln.de}
}
\note{
  Performance of the implemented methods has been evaluated and compared. The min-test performs very
  conservative if the means in the marginal treatment groups are close
  for a combination.
}
\seealso{\code{\link{bifactorial}}, \code{\link{carpet}}, \code{\link{cube}}, \code{\link{avetest}}, \code{\link{maxtest}},}
\examples{

#AML example from Huang et al. (2007) with data from
#Issa et al. (2004) and Petersdorf et al. (2007)
n<-c(10,31,17,100,50,50,101,50,50)
p<-c(0.00,0.45,0.65,0.30,0.71,0.70,0.59,0.64,0.75)
y<-list()
for(i in 1:9){
  y[[i]]<-0
  while((sum(y[[i]])!=round(n[i]*p[i]))||(length(y[[i]])==1)){
    y[[i]]<-rbinom(n[i],1,p[i])
  }
}

aml<-carpet(data=y,D=c(2,2))
mintest(aml,test="ztest",nboot=25000)
margint(aml,test="ztest",nboot=25000)




}
\keyword{multivariate}
\keyword{htest}
\keyword{iteration}
\keyword{design}
