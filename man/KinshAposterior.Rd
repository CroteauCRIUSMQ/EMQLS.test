\name{KinshAposterior}
\alias{KinshAposterior}
\title{
Kinship matrix a posterior
}
\description{
Function computes the a posterior kinship matrix.
}
\usage{
KinshAposterior(vec.prob, fam.id, sujet.id1, sujet.id2, data)
}
\arguments{
  \item{vec.prob}{Vector of the probability names in the database}
  \item{fam.id}{The family identifying}
  \item{sujet.id1}{The variable name of the first individual}
  \item{sujet.id2}{The variable name of the second individual}
  \item{data}{Data frame of the database.}
}
\details{
The function computes the a posterior kinship matrices according to the probabilities of identity by descent (IBD). For details, see \emph{T.Thornton et al., (2007)}.
}
\value{
Return a list contains the family kinship matrices.
}
\references{
Thornton T and McPeek MS: Case-control association testing with related individuals: a more powerful quasi-likelihood score test. \emph{Am J Hum Genet} 2007, 81:321, 337.
}
\seealso{
to See Also as \code{\link{kinshAprior}}
}
\examples{
 data(DataIBD)
    vec.prob=c("P0","P1","P2")
    fam.id="FAMILY"
    sujet.id1="ID1"
    sujet.id2="ID2"  
    lst_mkinshp1=KinshAposterior(vec.prob,fam.id,sujet.id1,sujet.id2,DataIBD)
    }