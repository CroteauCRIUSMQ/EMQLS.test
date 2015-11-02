flstkinsh <-
function(fam.id,sujet.id,dad.id,moth.id,sex,data,vec.prob){
          # this function permits to build the liste of kinship matrix of the family
          # fam.id: family identity
          # sujet.id: sujet identity
          # dad.id: dad identity
          # moth.id: mother identity
          # sex: Mal
          # data: database
          # vec.prob: probability vector 
          if(missing(vec.prob)==TRUE){
          lstMksp<-kinshAprior(fam.id,sujet.id,dad.id,moth.id,sex,data)
          }else{
                lsk1<-KinshAposterior(vec.prob,fam.id,sujet.id,data)
                lsk2<-kinshAprior(fam.id,sujet.id,dad.id,moth.id,sex,data)
                lstMksp<-list(lstkspAp=lsk2,lstkspApo=lsk1)
                }
          return(lstMksp)
          }
