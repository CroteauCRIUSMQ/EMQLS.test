\name{Mqls.test2}
\alias{Mqls.test2}
\title{
Extension MQLS test
}
\description{
This function computes the Extension MQLS test proposed by \emph{M. Nguile-Makao et al., (under review)}.
}
\usage{
Mqls.test2(data, lstk, fam.id, fl, k1, k)
}
\arguments{
  \item{data}{Data frame.}
  \item{lstk}{List contains the family kinship matrix.}
  \item{fam.id}{The family identifying}
  \item{fl}{Formula}
  \item{k1}{The nuisance parameter which represents the allelic frequency of the SNP associated to the endophenotype. For more detail see the  \emph{N. Nguile-Makao et al., (under review)}.}
  \item{k}{The nuisance parameter which represents the prevalence of the interaction between phenotype, endophenotype and SNP or the interaction between phenotype and SNP or endophenotype and SNP. By default the function computes the test statistic using the Davies's method. Remark: This is the SNP associated to the endophenotype.}
}
\details{
The function computes the extension of the MQlS proposed by \emph{M. Nguile-Makao et al (under review)}
}
\value{
Return a list containing

\item{Rstt }{The result matrix containing the statistic of the test, the number of freedom degree, the p-value computed with the parameter k and the corrected p-value computed by the Davies's method.}
\item{Rstd }{The result matrix containing the data description, the number of stratums, the number of freedom degree, the empiric allelic frequency, the estimated of the allelic frequency computed by BLUE method and the allelic frequency variance.}
\item{sg }{The number of kinship matrix that are singulars.}
\item{fMqls }{The function of the test statistic which depends of the nuisance parameter k.}
}
\references{
Moliere Nguile-Makao, Alexandre Bureau. Extension of the MQLS test: Application at Joint modelling of the complex diseases and endophenotypes to characterize genes and their interactions \emph{under review}.
	
Robert B. Davies. Hypothesis testing when a nuisance parameter is present only under the alternative \emph{Biometrika} (1987), 74, 1, pp.33-43.
}
\seealso{
To See Also as \code{\link{Mqls.test2}}
}
\examples{
data(DataGPh)
   lst_mkinshp2=kinshAprior("fam.id","subject.id","father.id","mother.id","sex",DataGPh)
   fl=formula(X2~Y2:X1)
   lR<-Mqls.test2(DataGPh,lst_mkinshp2,"fam.id",fl,0.3,0.38)
   lR
}

