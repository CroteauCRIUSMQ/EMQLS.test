\name{DataIBD}
\alias{DataIBD}
\docType{data}
\title{
Data frame of identity by descent (IBD)
}
\description{
Data frame contains the IBD probabilities of the target SNP. These database come from Alexandre works see \emph{A. Bureau et al., (2014)}.
}
\usage{data(DataIBD)}
\format{
  Data frame contains 13600 observations and 7 variables.
  \describe{
    \item{\code{FAMILY}}{ a numeric vector represents the family identifying.}
    \item{\code{ID1}}{ a numeric vector represents the individual 1 identifying.}
    \item{\code{ID2}}{ a numeric vector represents the individual 2 identifying.}
    \item{\code{MARKER}}{ a factor with levels \code{SNP_1} represents the name of the genetic marker.}
    \item{\code{P0}}{ a numeric vector represents the probability that the individual ID 1 and 2 share zero allele coming from of a common ancestor.}
    \item{\code{P1}}{ a numeric vector represents the probability that the individual ID 1 and 2 share one allele coming from of a common ancestor.}
    \item{\code{P2}}{ a numeric vector represents the probability that the individual ID 1 and 2 share two alleles coming from of a common ancestor.}
  }
}
\details{
The data was generated with the Merlin,QTDT format:
Add a URL: \url{http://csg.sph.umich.edu/abecasis/QTDT/docs/ibd.html} and convert them into R objects for the details see \emph{A.Bureau et al., (2014)}
}
\references{
Alexandre Bureau, Jordie Croteau, Yvon C. Chagnon, Marc-Andre Roy and Michel Maziade. Extension of the generalized disquilibrium test to polytomous phenotypes and two-locus models.
}
\examples{
data(DataIBD)
head(DataIBD)
}
