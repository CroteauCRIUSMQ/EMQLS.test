fZtild1 <-
function(datf,fl){
	    # the retrospective variable 
	    # typ: vector contains the type of the variable. Forma : the number 1 first case, 2 the second case, 3 the third case.
        # and the value of k1
  	    # this function returns the function that depends to parameter k2
        # step 1 built functions
        
        fpar1<-Vectorize(function(s){return(-s/(1-s))},'s')
        
        # step 1 builts the function that transforms the data in function according to k
        ftr.ztid<-function(vztild){
                  fvztild<-Vectorize(function(u){
                           return(ifelse(vztild==1,1,0)+ifelse(vztild==0,fpar1(u),0))
                           },'u')
                  return(fvztild)
                  }
        
         # Step 2 datamining
         vvar=all.vars(fl)[-1];repo<-all.vars(fl)[1]
         Ztildf<-as.vector(model.matrix(formula(paste("~",attr(terms.formula(fl),"term.labels"))),data=as.data.frame(apply(datf[vvar],2,function(u)ifelse(is.na(u)==TRUE,99,u))))[,2])
         zx2<-as.vector(datf[[repo]])
         X2<-as.vector(na.omit(zx2))
         # step 3 builts the Ztild
           if(dim(na.omit(datf))[1]==dim(datf)[1]){
                 fztild_m<-Vectorize(function(k){return(0)},'k');
                 fztild_n<-ftr.ztid(Ztildf[is.na(zx2)!=TRUE])
                 }else{fztild_m<-ftr.ztid(Ztildf[is.na(zx2)==TRUE])
                       fztild_n<-ftr.ztid(Ztildf[is.na(zx2)!=TRUE])
                      }
                 
         return(list(fztild_m=fztild_m,fztild_n=fztild_n,X2=X2))
  	   }
