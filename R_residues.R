library("R.matlab")

# The R script needs to be in the same folder as the .mat files.
setwd(getwd())

pDCont=readMat("prdDataControls.mat")$pD[,,1]
pDExp=readMat("prdDataExposed.mat")$pD[,,1]


for (name in names(pDCont)){
  temp=as.data.frame(pDCont[name])
  temp$diff=(temp[,2]-temp[,1])/temp[,2]*100        # (real values - prediction) / real values *100
  plot(temp$diff~temp[,3], main=name, type="l")
  abline(h=0, col="red", lty=2)
}


for (name in names(pDExp)){
  temp=as.data.frame(pDExp[name])
  temp$diff=(temp[,2]-temp[,1])/temp[,2]*100
  plot(temp$diff~temp[,3], main=name, type="l")
  abline(h=0, col="red", lty=2)
}
