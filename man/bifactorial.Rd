\name{bifactorial}
\alias{bifactorial-package}
\alias{bifactorial}
\docType{package}
\title{
  General information on the package
}
\description{
  Factorial clinical trial designs can be used to test for the efficacy of combination drugs with two or more components, where
  inference on the question if a combination therapy is more efficacious
  than both of its components is based on the min-test proposed by Laska and Meisner (1989). This is
  due to regulatoric demands requiring a contribution of all
  compounds in a combination drug. The AVE- and MAX-approaches
  proposed by Hung, Chi and Lipicky (1993) test for the existence of any desirable
  combination.
  
  Bootstrap-based methods are implemented as well as classical approaches available from
  literature to obtain p-values and confidence intervals in
  such designs. For the min-test, analytical methods use a normality and homoscedasticity
  assumption on the data (Hung, Chi and Lipicky, 1993 and Hung, 2000). Critical
  values needed for determination of confidence intervals are calculated using quantiles
  of the multivariate t-distribution (Bretz, Genz and Hothorn
  2001). These methods fail when handling with data that are skewed or
  heteroscedastic over the treatment groups. Furthermore, no analytical
  approach is available for the trifactorial case and the AVE- and
  MAX-tests on binary data.
  In the bootstrap approach, only the empirical distribution of the data is
  used and thus the results are valid for any distributional shape,
  provided that sufficiently large samples are available. Less
  analytical framework is needed to handle with the distributional
  properties of the tests. Further information on resampling-based
  methods and theoretical backgrounds are given in Westfall and Young
  (1993).

  Anyway, the problem of the extremely decreasing power
  for small values of the so-called nuisance parameters indicating
  the response differences between the marginal treatment groups cannot
  be resolved by the bootstrap approach. Any
  algorithm based on estimates for the nuisance parameters other than
  the assumption that they are infinite will exceed the given significance level (Snapinn, 1987).
    
  The package contains the generic functions \code{mintest} and \code{margint} to test for mean differences of given numeric data
  vectors and differences in event rates for binary data
  applications. Method dispatch is available for objects of class
  \code{carpet} or \code{cube}, which will lead to min-test
  results on a bi- or trifactorial design and corresponding confidence
  intervals comparing combination treatments with their respective
  component therapies. Implementations for global tests are also
  available by the generic functions \code{avetest} and \code{maxtest}.
}
\author{Peter Frommolt, University of Cologne \email{peter.frommolt@uni-koeln.de}\cr
  \url{http://www.medizin.uni-koeln.de/kai/imsie/homepages/peter.frommolt/}
}
\references{
  Bretz F, Genz A, Hothorn LA (2001): On the numerical availability of
  multiple comparison procedures. Biometrical Journal 43/5, pp. 645-656
  
  Hellmich M, Lehmacher W (2005): Closure procedures for monotone bi-factorial
  dose-response designs. Biometrics 61, pp. 269-276
  
  Hung HMJ, Chi GYH, Lipicky RJ (1993): Testing for the existence of a
  desirable dose combination. Biometrics 49, pp. 85-94

  Hung HMJ, Wang SJ (1997): Large-sample tests for binary outcomes in
  fixed-dose combination drug studies. Biometrics 53, pp. 498-503
  
  Hung HMJ (2000): Evaluation of a combination drug with multiple doses in
  unbalanced factorial design clinical trials. Statistics in Medicine
  19, pp. 2079-2087

  Laska EM, Meisner MJ (1989): Testing whether an identified treatment is
  best. Biometrics 45, pp. 1139-1151

  Snapinn SM (1987): Evaluating the efficacy of a combination
  therapy. Statistics in Medicine 6, pp. 657-665
  
  Westfall PH, Young SS (1993): Resampling-based multiple testing. John
  Wiley & Sons, Inc., New York
}
\keyword{package}
\keyword{design}
\keyword{htest}
\keyword{multivariate}
\keyword{iteration}
