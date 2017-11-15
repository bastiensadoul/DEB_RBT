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
  } else {
    list(c(-ks*y[1]/cs,
           0)
    )
  }
  
}





############################ THE MODEL

# Test other functions

exp_decrease=function(t,y, param){
  
  ks <- param[1]
  cs <- param[2]
  
  if (t>tmin && t<tmax){
    list(dy = 0)
  } else {
    if (identical_recovery_time == "TRUE"){list(dy = -y*1/ks)     # used to be -y*ks when "best fits manually" was created
    } else {list(dy = -y/Fpert * 1/ks)} 
  }
  
}


decreasing_logistic=function(t,y, param){
  
  ks <- param[1]
#  cs <- param[2]
  
  if (t>tmin && t<tmax){
    list(dy = 0)
  } else {
    list(dy = -ks*y*(1-y))  # can divide the second y by cs to add another parameter
  }
  
}


linearmod=function(t,y, param){
  
  ks <- param[1]
  #  cs <- param[2]
  
  if (t>tmin && t<tmax){
    list(dy = 0)
  } else if (y>0) {
      if (identical_recovery_time == "TRUE"){list(dy = -Fpert/ks)
        } else {list(dy = -1/ks)}
  } else {list(dy = 0)}
  
}



