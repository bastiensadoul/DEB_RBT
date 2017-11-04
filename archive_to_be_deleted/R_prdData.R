library("R.matlab")

# You need to change "dir", to make that script work (address of the folder of interest)
dir="~/GitHub/DEB_Oncorhynchus_mykiss/DEB_RBT"
setwd(dir)


res=readMat("prdData.mat")$prdData[,,1]



for (name in names(res)){
  temp=as.data.frame(res[name])
  temp$diff=(temp[,3]-temp[,2])/temp[,3]*100        # (real values - prediction) / real values *100
  plot(temp$diff~temp[,1], main=name, type="l")
  abline(h=0, col="red", lty=2)
}
