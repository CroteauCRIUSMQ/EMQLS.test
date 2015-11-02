GenoCod_Gen <-
function(data,MA,genotyp,Method=c("recessive","dominant","codominant"),var.nam=NULL,na.rm=FALSE){
          # Method: gives the method to use 
          # data: database 
          # na.rm: treatment (TRUE) or not of the missing values (False) 
          # MA: minor allele
          # genotyp
          # this function permits to count the number of the minor allele by individual
          
          Method=match.arg(Method)
          switch(Method,recessive=GenoCod_Gr(data,MA,genotyp,var.nam,na.rm),
                        dominant=GenoCod_Gd(data,MA,genotyp,var.nam,na.rm),
                        codominant=GenoCod_Gcod(data,MA,genotyp,var.nam,na.rm))
          }
