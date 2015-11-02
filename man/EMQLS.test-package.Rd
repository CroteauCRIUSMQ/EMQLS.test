\name{EMQLS.test-package}
\alias{EMQLS.test-package}
\alias{EMQLS.test}
\docType{package}
\title{
Extension of Maximum Quasi-likelihood Score test
}
\description{
Package implements the Maximum Quasi-likelihood Score test proposed by \emph{T. Thornton et al., (2004)} and an extension of the test proposed by \emph{M. Nguile-Maka et al., (under review)}.

}
\details{
\tabular{ll}{
Package: \tab EMQLS.test\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-04-02\cr
License: \tab GPL-2\cr
}
The package \pkg{EMQLS.test} implements the Maximum Quasi-likelihood Score test proposed by \emph{T. Thornton et al., (2004)} and an extension of the test proposed by \emph{M. Nguile-Maka et al.,(under review)}. EMQLS permits to detect the association between SNPs and the complex phenotypes in case-control samples where individuals are related in the stratum (Families). The extension proposed permits to include the endophenotype and her SNP associated, in the aim to increase the test power. The package also proposes other functions like the function to compute the posterior kinship matrices \code{KinshAposterior}, minor allelic frequency via the best linear unbiased estimator (BLUE) (\code{Beta_est}) EMQLS test (\code{MQLS.test1},\code{MQLS.test2}). Package uses the kinship function of the package \pkg{kinship2} to estimate the a posterior kinship matrix and the function \code{pseudoinverse} of the package \pkg{corpcor} to estimate the pseudo-inverse the kinship matrix when it has the problem of the singularity matrix.
}
\author{
Moliere Nguile Makao & Alexandre Bureau.

Maintainer: Moliere Nguile Makao <moliere.nguile@gmail.com>
}
\references{
Timothy Thornton and Mary sara McPeek Case-Control Association Testing with Related Individuals: A More Powerful Quasi-likelihood Score Test \emph{American journal of human genetics} vol 81 August 2007

Alexandre Bureau, Jordie Croteau, Yvon C. Chagnon, Marc-Andre Roy and Michel Maziade. Exension of the generalized disequilibrium test to polytomous phenotypes and two-locus models \emph{Genetics} DOI: 10.3389/fgene.2014.00258.

Moliere Nguile-Makao, Alexandre Bureau. Extension of the MQLS test: Application at Joint modelling of the complex diseases and endophenotypes to characterize genes and their interactions \emph{under review}
}
\seealso{
\code{\link{Beta_est}}, \code{\link{Deco.kinf}}, \code{\link{GenoCod_Gen}}, \code{\link{KinshAposterior}}, \code{\link{kinshAprior}}, \code{\link{Mqls.test1}}, \code{\link{Mqls.test2}}
}

