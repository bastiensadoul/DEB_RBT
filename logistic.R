logistic=function(x, param){
  k=param$k
  xo=param$xo
  hipM=param$hipM                                       # percentage of increase of pM at the first dpf
  pM=param$pM
  return((hipM*pM/100)-pM/(1+exp(-k*(x-xo))))
}

temp=seq(0,500, length.out=500)

param=data.frame(k=0.02
                 ,xo=150
                 ,hipM=200
                 ,pM=173                               # this should be the value of pM in the control
                 )

pM_pred=logistic(temp, param)

plot(pM_pred~temp, type="l", col="red")
abline(h=param$pM, lty=2)
