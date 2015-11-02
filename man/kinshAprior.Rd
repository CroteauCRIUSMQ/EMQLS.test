\name{kinshAprior}
\alias{kinshAprior}
\title{
Kinship matrix a prior
}
\description{
The function computes and creates a list containing of the a prior kinship matrices of the families.
}
\usage{
kinshAprior(fam.id, sujet.id, dad.id, moth.id, sex, data)
}
\arguments{
  \item{fam.id}{Family identifying}
  \item{sujet.id}{Individual  identifying within the family}
  \item{dad.id}{Name of father variable.}
  \item{moth.id}{Name of mother variable.}
  \item{sex}{Name of sex variable.}
  \item{data}{Data frame.}
}
\details{
The function uses the function \code{\link{kinship}} of the package \pkg{kinship2} to assess a prior family kinship matrices.
}
\value{
Return a list of the a prior family kinship matrices.
}
\seealso{
to See Also as \code{\link{KinshAposterior}}
}
\examples{
data(DataGPh)
lst_mkinshp2=kinshAprior("fam.id","subject.id","father.id","mother.id","sex",DataGPh)
}

