\name{Mqls.test1}
\alias{Mqls.test1}
\title{
Computed the MQLS test
}
\description{
This function computes the MQLS proposed by Thormton. see the details \emph{T. Thornton et al. 2004}.
}
\usage{
Mqls.test1(data, lstk, fam.id, fl, k)
}
\arguments{
  \item{data}{Data frame.}
  \item{lstk}{List of the family kinship matrices}
  \item{fam.id}{Family identifying}
  \item{fl}{The formula}
  \item{k}{The disease prevalence. By default, it's equal to zero. A second statistic is evaluated while optimizing the parameter k by Davies method. see \emph{D. Davies et al. (1987)}
}
}
\details{
The fonction implantes the MQLS published by \emph{T. Thornton et al (2004)}. The function uses the Davies method \emph{D. Davies et al. (1987)} to evaluate a second statistic while optimizing the parameter k. A Software in C++ was also developed for the test. see \url{https://galton.uchicago.edu/~mcpeek/software/MQLS/index.html}.
}
\value{
Return a list containing

\item{Rstt }{The result matrix containing the statistic of the test, the number of freedom degree, the p-value computed with the parameter k and the corrected p-value computed by the Davies's method.}
\item{Rstd }{The result matrix containing the data description, the number of stratums, the number of freedom degree, the empiric allelic frequency, the estimated of the allelic frequency computed by BLUE method and the allelic frequency variance.}
\item{sg }{The number of kinship matrix that are singulars.}
\item{fMqls }{The function of the test statistic which depends of the nuisance parameter k.}
}
\references{
Timothy Thornton and Mary sara McPeek Case-Control Association Testing with Related Individuals: A More Powerful Quasi-likelihood Score Test \emph{American journal of human genetics} vol 81 August 2007.
	
	
Robert B. Davies. Hypothesis testing when a nuisance parameter is present only under the alternative \emph{Biometrika} (1987), 74, 1, pp.33-43.
}

\seealso{
To See Also as \code{\link{Mqls.test2}}
}
\examples{
	data(DataGPh)
    fl=formula(X2~Y2)
    lst_mkinshp2=kinshAprior("fam.id","subject.id","father.id","mother.id","sex",DataGPh)
    lR<-Mqls.test1(DataGPh,lst_mkinshp2,"fam.id",fl,0.008)
    lR
}
