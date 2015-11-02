kinshAprior <-
function(fam.id,sujet.id,dad.id,moth.id,sex,data){
               # fam.id: family identity
               # sujet.id: sujet identity
               # dad.id: dad variable 
               # moth.id: mother variable
               # sex: mal variable 
                 vec.fam<-unique(data[[fam.id]])
               # Creation de la fonction kinship
                 fksh<-function(u){dat<-data[data[fam.id]==u,]
                                   if(dim(dat)[1]==1){stop("Familly having one individual")
                                      }else{
                                            return(2*kinship(dat[[sujet.id]],dat[[dad.id]],dat[[moth.id]],dat[[sex]]))
                                            }
                                   }
               # creation de la liste de kinship 
                 lskin<-lapply(vec.fam,fksh)
                 names(lskin)<-vec.fam
                 return(lskin)
                 }
