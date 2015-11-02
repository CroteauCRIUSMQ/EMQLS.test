\name{Beta_est}
\alias{Beta_est}
\title{
Computes the allelic frequency
}
\description{
Computes the allelic frequency using the Best Linear Unbiased Estimator (BLUE) and the kinship matrices
}
\usage{
Beta_est(lstfki)
}
\arguments{
  \item{lstfki}{
MQLS object computed by the function \code{\link{Deco.kinf}}.
}
}
\details{
The function uses the Best Linear Unbiased Estimator (BLUE) and the family kinship matrices to compute the allelic frequency. For the details, see \emph{T. Thornton et al.,(2007).}
}
\value{
Return the allelic frequency.
}
\references{
Thornton T and McPeek MS: Case-control association testing with related individuals: a more powerful quasi-likelihood score test. \emph{Am J Hum Genet} 2007, 81:321,337.

Bourgain C, Hoffjan S, Nicolae R, Newman D, Steiner L, Walker K, Reynolds R, Ober C and McPeek MS:Novel case-control test in a founder population identifies P-selectin as an atopy-susceptibility locus. \emph{Am J Hum Genet} 2003, 73:612,626.
}
\seealso{
\code{\link{Deco.kinf}}
}
