KinshAposterior <-
function(vec.prob,fam.id,sujet.id1,sujet.id2,data){
               # sujet.id1 it must that in database the sujet.id is eqal to sujet.id2
               # vec.prob=c(P0,P1,P2) probability vector of the copi number.
               # the probabilty to sharing i allele(s) with another family mumber.
                 p<-length(vec.prob)
                 
               # la petite fonction qui traite les données maquantes
                 ftm<-function(vec){
                                     vep<-c(0,0.5,1)
                                     vec1<-as.vector(vec)
                                    if(length(vec1[is.na(vec1)!=TRUE])>1){
                                        if(length(vec1[is.na(vec1)!=TRUE])==2){
                                           vec2<-vec1
                                           vec2[which(is.na(vec2)==TRUE)]<-(1-sum(vec1[is.na(vec1)!=TRUE]))
                                                                            }else{vec2<-vec1}
                                                                        }else{
                                                                            if(sum(vec1[is.na(vec1)!=TRUE])==1){
                                                                            vec2<-vec1
                                                                            vec2[which(is.na(vec2)==TRUE)]<-0
                                                                            }else{vec2<-vec1}
                                                                        }
                                    return(vec2%*%t(t(vep)))}
                                             
                 ksh.p<-apply(as.matrix(data[vec.prob]),1,ftm)
                 dat1<-data[-c(which(names(data)%in%vec.prob))]
                 datalf<-data.frame(dat1,ksh.p)
               # creation de la matrice de kinship 
                 vec.fam<-unique(datalf[[fam.id]]) 
                 lstfam<-lapply(vec.fam,function(u)datalf[datalf[fam.id]==u,]) 
                 names(lstfam)<-vec.fam
                
               # la fonction qui transforme une liste en matrice
                 listM<-function(vec){uu<-max(sapply(vec,length))
                       return(t(sapply(vec,function(u) c(u,rep(NA,uu-length(u))))))} 
                
               # la liste de matrice de kinship appriorie            
                 lstk<-lapply(vec.fam,function(u){mat<-listM(lapply(unique(lstfam[[as.character(u)]][[sujet.id2]]),
                             function(v) t(lstfam[[as.character(u)]][lstfam[[as.character(u)]][sujet.id1]==v,]["ksh.p"])));
                             mat[upper.tri(mat)]<-mat[lower.tri(mat)]
                             return(mat)}) 
                 names(lstk)<-vec.fam                            
                 return(lstk)
                 }
