####################################################################################################################################################
####     Developped by Sadoul and co-workers for the software R. See :
####     "On the use of a simple physical system analogy to study robustness features in animal sciences"
####     Plos One 2015
####################################################################################################################################################



############################ THE MODEL

# This function gives for a given time t, for a given state y of the system (vector containing the speed and the acceleration at time t) 
#        and of the parameters param (vector of ks and cs) 
# ks and cs are the properties of the spring and damper respectively

spring_damper_model=function(t,y, param){
  
  ks <- param[1]
  cs <- param[2]
  if (length(param)==3){
    perc <- param[3]
  } else perc=1
  
  if (t>tmin && t<tmax){
    list(c(y[2],
           1/(perc*Fpert)*(perc*Fpert-ks*y[1]-cs*y[2])
    ))
  }
  else {
    list(c(-ks*y[1]/cs,
           0)
    )
  }
  
}