Beta_est <-
function(lstfki){
	      # This function estimated the frequency of the minor allele
  	      # lstfki: the object of the function Deco.kinf
          
          # Calculation of the two compoments
          Phi_n<-lstfki[["phi_n"]];Phi_m<-lstfki[["phi_m"]]
          Phi_ninv<-lstfki[["phi_ninv"]];X2<-lstfki[["X2"]]
          pa1<-0;kof<-0;
          idfam<-names(Phi_n)
  	      for(obs in idfam){
  	                        Xb2<-X2[[obs]]
  	                        Dp<-rep(1,length(Xb2))
  	                        # Estimation
  	                        pa1<-pa1+t(Dp)%*%Phi_ninv[[obs]]%*%Dp  # Dp: the difference between the expectation and Beta
  	                        kof<-kof+t(Dp)%*%Phi_ninv[[obs]]%*%Xb2  # Phi_ninv: The inverse of the kinship matrix
  	                        }
         P.est<-as.numeric(kof/pa1)
         return(P.est)
         }
