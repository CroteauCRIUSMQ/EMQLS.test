GenoCod_Gcod <-
function(data,MA,genotyp,var.nam=NULL,na.rm=FALSE){
          # Mode codominant 
          # data: database 
          # na.rm: treatment (TRUE) or not of the missing values (False) 
          # MA: minor allele
          # genotyp
          # this function permits to count the number of the minor allele by individual
          # creation of observation vector
            obs<-c(1:dim(data)[1])
            data1<-data.frame(obs,data)
          # verification of arguments
          if(is.factor(data1[[genotyp]])!=TRUE){stop("It's not a character")}
          else{
               # replacement empty spaces by the NA 
                 data1[genotyp][data1[genotyp]==""]<-NA
                 if(na.rm==TRUE){data1<-data1
                 }else{data1<-data1[is.na(data1[genotyp])!=TRUE,]}
                      # recover the ID
                        Id1<-data1["obs"] 
                      # suppression of the character that aren't the letters 
                        uu1<-lapply(apply(data1[genotyp],2,strsplit,"")[[1]],function(vv)vv[vv%in%letters|vv%in%LETTERS])
                            
                      # Transformation the liste at the matrix 
                        listM<-function(vec){uu<-max(sapply(vec,length))
                               return(t(sapply(vec,function(u) c(u,rep(NA,uu-length(u))))))} 
                                        
                      # split the matrix
                        mat1<-listM(uu1) 
                        n<-dim(mat1)[2]
                        n1<-n/2;vec1<-c(1:n1);
                        n2<-n1+1;vec2<-c(n2:n)
                        mtge1<-mat1[,vec1]
                        mtge2<-mat1[,vec2]
               
                      # The function that counts the alleles
                        ff<-function(zz){
                        	 if(sum(is.na(zz))==n1){zz1<-NA
                              }else{
                                    if(0<sum(zz%in%MA) & sum(zz%in%MA)<n1+1){zz1<-1}else{zz1<-0}
                                    }
                             return(zz1)
                            }
                    # les derniers calculs de la fonction
                      Allel_p<-apply(mtge1,1,ff)
                      Allel_m<-apply(mtge2,1,ff)
                    # Calculation of the allelic frequency
                       alf<-(Allel_p+Allel_m);
                      p<-sum(alf,na.rm=TRUE)/(2*(length(Allel_m[!is.na(Allel_m)])))
                      Allel<-(Allel_p+Allel_m)/2 
                      datR<-as.data.frame(cbind(Id1,Allel))
                      if(is.null(var.nam)==TRUE){
                      names(datR)<-c("obs","X1")
                      }else{names(datR)<-c("obs",var.nam)}
                      #data1[genotyp]<-NULL
                      basR<-merge(data1,datR,by=c("obs"))
                      basR["obs"]<-NULL
                    # sortite 
                      return(list(datR=basR,MA.frq=p))
                      }
                }
