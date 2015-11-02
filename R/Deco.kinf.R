Deco.kinf <-
function(dataf,matkaf,fam.id,fl){
  	    # The function to define the phi matrix
  	    # Note when it hasn't the missing values, the phi is reduced to kinship
  	    # The same function also permits to compute the kinship inverse according the several cases                                           	
  	    # datf: the familly database
  	    # matkf: The familly kinship matrix
        fam.obs<-unique(dataf[[fam.id]])
        
        # initialisation
        phi_n=list();phi_ninv=list();
        phi_m=list();X2=list();sg=0
        vvar<-all.vars(fl);lst_datf<-list()
        
        #starting the program
        for(obs in fam.obs){
            datf<-dataf[dataf[[fam.id]]==obs,];matkf<-matkaf[[as.character(obs)]]
            lst_datf[[as.character(obs)]]<-datf
        if(is.null(dim(matkf))==TRUE){stop("Error: the number of kinship matrix is upper to the number of family")
        }else{
  	      Xb2<-datf[is.na(datf[[vvar[1]]])!=TRUE,][[vvar[1]]]    # on supprime les valeurs manquantes
  	      n<-dim(datf[is.na(datf[[vvar[1]]])!=TRUE,])[1];        # nombre d'individus n'ayant pas des genotypes manquants
  	      m1<-as.numeric(dim(datf)[1]);                          # nombre total des membres de la famille
  	      m<-as.numeric(m1-n);v1<-n+1
  	      # If it hasn't the missing genotypes  
          if(m==0){phi_n0<-matkf;phi_n[[as.character(obs)]]<-phi_n0;
                   phi_m[[as.character(obs)]]<-0
                   }else{
                   phi_n0<-matkf[1:n,][,1:n];phi_n[[as.character(obs)]]<-phi_n0;
                   phi_m[[as.character(obs)]]<-as.matrix(matkf[1:n,])[,v1:m1]}
  	      # calcul de la matrice inversible 
  	      if(is.character(try(solve(phi_n0),TRUE))==TRUE){phi_ninv[[as.character(obs)]]<-pseudoinverse(phi_n0);sg=sg+1
  	      	                        }else{phi_ninv[[as.character(obs)]]<-solve(phi_n0);sg=sg}
          X2[[as.character(obs)]]<-Xb2
        }
        }
  	      lstK<-list(phi_n=phi_n,phi_ninv=phi_ninv,phi_m=phi_m,X2=X2,lst_datf=lst_datf,sg=sg)
  	      return(lstK)
        }
