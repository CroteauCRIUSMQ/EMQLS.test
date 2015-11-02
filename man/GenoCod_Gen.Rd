\name{GenoCod_Gen}
\alias{GenoCod_Gen}
\title{
Transforms the genetic maker
}
\description{
Function transforms the genetic marker to recessive, dominant or codominant mode.
}
\usage{
GenoCod_Gen(data, MA, genotyp, Method, var.nam = NULL, na.rm = FALSE)
}
\arguments{
  \item{data}{Data frame of the database}
  \item{MA}{Minor allele name in the database}
  \item{genotyp}{Genotype name in the database}
  \item{Method}{The method to use\code{Method = c("recessive", "dominant", "codominant")}}
  \item{var.nam}{Name of the new variable, by default the name is \code{X2}}
  \item{na.rm}{TRUE, the function conserves the missing data and FALSE, it suppresses the individuals, which have the missing data.}
}
\details{
Function that transforms the genetic marker according to the mode ("recessive", "dominant", "co-dominant") and assess the allele minor frequency under Hardy Weinberg equilibrium of the target population.
}
\value{
Return a list containing
\item{datR}{Data frame with the new variable}
\item{MA.frq}{Minor allele frequency}
}
\note{
We use the documentation of the site below to code the genetic markers in "Dominant", "Recessive" and "codominant" mode.
}
\references{
\url{http://wims.unice.fr/wims/wims.cgi?module=U1/genetics/docpopulation.fr}
}
\examples{
# Example : Genotype data

  set.seed(123000)
  snp1 <- c("aaAA","AA:aA","Aa/AA","aa/AA","aA/aa","Aa/Aa","AA/Aa","AA/Aa","AA/Aa",NA)
  snp2 <- c("Aa/aA","Aa/AA","AA/aa","AA/aA","AA/Aa","Aa,AA","AA/AA","AA/aA","Aa/aa","")
  datD<-NULL
  for(u in 1:10){
              SNP1<-sample(snp1)
              SNP2<-sample(snp2)
              var1 <- runif(length(snp1))
              pheno<- rbinom(length(snp1),c(0,1),prob=c(0.65,0.35))
              pheno[sample(1:length(snp1))[3]]<-NA
              id<-c(1:length(snp1))
              id.fam<-rep(u,length(snp1))
              tabGen1<-data.frame(id.fam,id,SNP1,SNP2,var1,pheno)
              datD<-rbind(datD,tabGen1)
              }
  # Result
    head(datD)
    
  # Example of the using GenoCod_Gen
  # 'a' is the minor allele, method "recessive" the new name is X2
   datT1<-GenoCod_Gen(datD,"a","SNP2",Method=c("recessive"),var.nam="X2")
   head(datT1$datR)
}

