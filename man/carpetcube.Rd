\name{carpetcube}
\alias{carpet}
\alias{carpet-class}
\alias{carpet,list,numeric-method}
\alias{carpet,ANY,ANY-method}
\alias{show,carpet-method}
\alias{summary,carpet-method}
\alias{plot,carpet,missing-method}
\alias{cube}
\alias{cube-class}
\alias{cube,list,numeric-method}
\alias{cube,ANY,ANY-method}
\alias{show,cube-method}
\alias{summary,cube-method}
\alias{plot,cube,missing-method}
\title{Objects for handling with bi- and trifactorial trial data}
\description{
  Create objects representing bi- or trifactorial clinical trial designs.
}
\usage{
carpet(data,D,...)
cube(data,D,...)
}
\arguments{
  \item{data}{A list of numeric or binary data vectors from the
    trial. See the details below for the order in which the list is to be given.}
  \item{D}{An integer vector of length 2 for \code{carpet} objects and
    of length 3 for \code{cube} objects, specifying the number of doses of the components drugs in the trial.}
  \item{\dots}{Any further arguments.}
}
\details{
  The function \code{carpet} creates objects of class \code{carpet} from the
  specified \code{data} in the list that are used row-wise to fill up
  the 2-factorial treatment groups, i.e. in the order
  (0,0), (0,1),..., (0,D[2]), (1,0), ..., (1,D[2]), ..., (D[1],D[2]);
  resulting in a (D[1]+1)x(D[2]+1) data array.
  
  To represent trifactorial designs for the evaluation of a
  three-compound combination, an object of class \code{cube} can be
  created using the function \code{cube}. The data in the treatment
  groups are then filled up in the order (0,0,0), ..., (0,0,D[3])
  first, then (0,1,0), ..., (0,1,D[3]) and up to
  (0,D[2],0), ..., (0,D[2],D[3]). This is the order also for
  the values 0, ..., D[1] for the first component group, always
  taking the data succesively from the list elements of \code{data}. The
  result is a (D[1]+1)x(D[2]+1)x(D[3]+1) data array. Methods for
  multiple inference and global tests can
  be applied to \code{carpet} and \code{cube} objects.
}
\value{
  An object of class \code{carpet} or \code{cube}, respectively, with
  the following slots.
  \item{data}{The data list specified in the construction.}
  \item{D}{Vector of maximum doses specified in the construction.}
  \item{n}{Numeric vector of sample sizes in the respective groups.}
}
\references{
  Hung HMJ, Chi GYH, Lipicky RJ (1993): Testing for the existence of a
  desirable dose combination. Biometrics 49, pp. 85-94

  Hung HMJ (2000): Evaluation of a combination drug with multiple doses
  in unbalanced factorial design clinical trials. Statistics in
  Medicine 19, pp. 2079-2087
}
\author{Peter Frommolt, University of Cologne \email{peter.frommolt@uni-koeln.de}\cr
  \url{http://www.medizin.uni-koeln.de/kai/imsie/homepages/peter.frommolt/}
}
\seealso{\code{\link{bifactorial}}, \code{\link{mintest}},
  \code{\link{margint}}, \code{\link{avetest}}, \code{\link{maxtest}}}
\examples{
#Hypertension example from Hung (2000)
data(sidbp)
x<-split(sidbp$ynrmhom,sidbp$cb)
bifactorial<-carpet(data=x,D=c(2,3))
}
\keyword{design}
\keyword{classes}
