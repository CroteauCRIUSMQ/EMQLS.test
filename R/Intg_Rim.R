Intg_Rim <-
function(fct,a,b,ndv=2){
               pas=(b-a)/ndv
               vecsd=seq(a,b,pas)
               vecsp<-vecsd[-1]
               vecinf<-vecsd[-length(vecsd)]
               return(sum(abs(fct(vecinf)-fct(vecsp))))
               }
