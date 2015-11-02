\name{DataGPh}
\alias{DataGPh}
\docType{data}
\title{
Simulated genetic data
}
\description{
Data frame that contains genetic variants with two alleles such as single nucleotide polymorphisms (SNPs) at two independent loci and Two phenotypes \eqn{Y1} and \eqn{Y2}.
}
\usage{data(DataGPh)}
\format{
  Data frame contains 100 families of 16 members and 11 variables.
  \describe{
    \item{\code{fam.id}}{  a numeric vector represents the family identifying.}
    \item{\code{subject.id}}{  a numeric vector represents the individual identifying within the family.}
    \item{\code{father.id}}{  a numeric vector represents the father identifying.}
    \item{\code{mother.id}}{  a numeric vector represents the mother variable.}
    \item{\code{sex}}{  a numeric vector represents the sex variable.}
    \item{\code{Y1}}{  a numeric vector of the endophenotype variable}
    \item{\code{Y2}}{  a numeric vector of the phenotype variable}
    \item{\code{X1}}{  a numeric vector of first loci variable}
    \item{\code{X2}}{  a numeric vector of the second loci variable}
    \item{\code{V14}}{  a numeric vector the allele 1 for X2}
    \item{\code{V15}}{  a numeric vector the allele 2 for X2}
  }
}
\details{
Data frame comes from the works of \emph{A.Bureau et al., (2014)}. The database contains 100 families of 16 members and 10 variables. It contains genetic variants with two alleles such as single nucleotide polymorphisms (SNPs). The genotypes founders were sampled under Hardy-Weinberg equilibrium using risk allele frequencies (RAFs) of 0.1 at locus 1 and 0.3 at locus 2. The transmission of alleles to their descendants was then simulated following the rules of mendelien inheritance. Two independent loci and two phenotypes \emph{Y1} and \emph{Y2}. For the details, see the simulation section of the article of  \emph{A.Bureau et al., (2014)}.
}
\references{
Alexandre Bureau, Jordie Croteau, Yvon C. Chagnon, Marc-Andre Roy and Michel Maziade. Extension of the generalized disequilibrium test to polytomous phenotypes and two-locus models. DOI 10.3389/fgene.2014.00258.
}
\examples{
data(DataGPh)
head(DataGPh)
}

