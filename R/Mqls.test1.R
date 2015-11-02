Mqls.test1 <-
function(data,lstk,fam.id,fl,k){
	# The Mqls test
	# data : the database 
	# lstk : the liste contain the kinship matrix of familly 
	# fam.id : the familly id 
	# fl : the formula
	# typ : the model type
	# k1 : the allelic frequecy
    
	# step2 to compute the score function
      lstfki<-Deco.kinf(data,lstk,fam.id,fl)
              Phi_n<-lstfki$phi_n;
              Phi_m<-lstfki$phi_m
              Phi_ninv<-lstfki$phi_ninv;
              X2<-lstfki$X2
              fam.obs<-names(Phi_n);
              slt_datf<-lstfki$lst_datf
              sg<-lstfki$sg
              
              # The calculation of Beta and rho
              Beta<-Beta_est(lstfki)
              ro2<-as.numeric((1/2)*(Beta*(1-Beta)))
    # the function Mqls
      fMqls<-Vectorize(function(k){
           # intial terms
           Ast<-0;Usc<-0;delta<-0;ct<-0;
           
           # the curly
           for(obs in fam.obs){
               datf<-slt_datf[[obs]]
               # the kinship matrix
               phi_n<-Phi_n[[obs]];
               phi_m<-Phi_m[[obs]]
               phi_ninv<-Phi_ninv[[obs]]
               ZX2<-X2[[obs]]
               Df<-as.vector(rep(1,length(ZX2)));
               
               # the variable Ztild
               Ztild<-fZtild1(datf,fl)
               # Compoment Xb2=(X2-mu)
               Xb2<-as.vector(ZX2-Df*Beta)
               # Compoment alpha
               if(is.null(dim(phi_m)[1])==TRUE){alph<-Ztild$fztild_n(k)}else{alph<-Ztild$fztild_n(k)+phi_ninv%*%phi_m%*%Ztild$fztild_m(k)}
               # Compoment omega
               if(is.null(dim(phi_m)[1])==TRUE){wf<-phi_n%*%Ztild$fztild_n(k)}else{wf<-phi_n%*%Ztild$fztild_n(k)+phi_m%*%Ztild$fztild_m(k)}
               # Compoment Usc
               Usc<-Usc+t(Xb2)%*%alph
               # Compoment A
               Ast<-Ast+t(alph)%*%wf
               # composante delta
               delta<-delta+t(alph)%*%Df
               # Compoment ct
               ct<-ct+t(Df)%*%phi_ninv%*%Df
           }
           # calculation of the statistic test
           # Gama Compoment
             Gama<-(Ast-(delta%*%t(delta))/ct)
           # Mqls statistic
             Mqls<-(Usc**2)/(ro2*Gama)
             return(Mqls)
       },'k')
      
       # step3 optimisation
       # on suppose que notre intervalle sera toujour c(0,1)
         sqfct<-Vectorize(function(x){sqrt(fMqls(x))^(1/2)},'x')
       
       # The optimize the function Mqls
         M1<-optimize(fMqls,interval=c(0,0.99),lower=0,upper=0.99,maximum=TRUE)$objective
       
       # To compute the integrate of the function
         V<-Intg_Rim(sqfct,0,0.99)
       
       # the corrected p-value
         Pco1=(1-pchisq(M1,1))+V*exp(-(1/2)*M1)*(2*pi)^(-1/2)
       
       
       if(missing(k)==TRUE){vdev<-c(round(fMqls(0),digits=4),1,round((1-pchisq(fMqls(0),1)),digits=4),round(Pco1,digits=4))
       }else{vdev<-c(round(fMqls(k),digits=2),1,round((1-pchisq(fMqls(k),1)),digits=4),round(Pco1,digits=4))}
       
       # the resultats
       ver<-t(t(t(vdev)))
       colnames(ver)<-c("Mqls","df","Pr(>|z|)","Pr(>|z|)*")
       rownames(ver)<-labels(terms(fl))
       
       # tabls stat
       name.cov=all.vars(fl);
       datD1<-data[is.na(data[[name.cov[1]]])!=TRUE,];# calcul du beta empirique
       Beta.emp=mean(datD1[[name.cov[1]]]);# calcul du beta empirique
       vcst<-t(t(t(c(length(fam.obs),Beta.emp,Beta,ro2))))
       colnames(vcst)<-c("N.strata","Beta.emp","Beta.est","Var.est")
       rownames(vcst)<-all.vars(fl)[1]
       
       return(list(Rstt=ver,Rstd=vcst,sg=sg,fMqls=fMqls))
      }
