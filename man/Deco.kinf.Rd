\name{Deco.kinf}
\alias{Deco.kinf}
\title{
MQLS object
}
\description{
The function creates the MQLS object.
}
\usage{
Deco.kinf(dataf, matkaf, fam.id, fl)
}

\arguments{
  \item{dataf}{Data frame.}
  \item{matkaf}{List that contains the family kinship matrices.}
  \item{fam.id}{Family identifying.}
  \item{fl}{Formula.}
}
\details{
The function prepares the necessary arguments to compute the MQLS test. It creates a list of the kinship matrices and their inverses. The kinship matrices are assessed while taking into account the missing data on target SNP. In case where there has the missing data, the function returns two matrices, one for the complete data and the other for incomplete data. For more detail see \emph{Nguile-Makao et al., (2015)}.
}
\value{
Return a list containing components,
  \item{phi_n}{List of the kinship matrices for n individual(s) not having the missing data on the target SNP.}
  \item{phi_ninv}{List of the inverse or pseudo-inverse kinship matrices for n individual(s).}
  \item{phi_m}{List of the kinship matrices for m individual(s) having the missing data on the target SNP.} 
  \item{label target marker}{List of the target genetic marker by family.}
  \item{lst_datf}{List of the Database by family.}
  \item{sg}{Number of the kinship matrices that are not invertible.}
}
\references{
Moliere Nguile-Makao, Alexandre Bureau (2015),Extension of the MQLS test: Application at joint modelling of complex diseases and endophenotypes to characterize genes effects and their interactions \emph{under review}.
}

\seealso{
To See Also as \code{\link{kinshAprior}}
}
\examples{
data(DataGPh)
lst_mkinshp2=kinshAprior("fam.id","subject.id","father.id","mother.id","sex",DataGPh)
fl=X2~X1:Y2;
fam.id="fam.id"
lst_fam21<-Deco.kinf(DataGPh,lst_mkinshp2,fam.id,fl)
}
