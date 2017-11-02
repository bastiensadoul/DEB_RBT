rm(list=ls())
cat("\014")  # To clear the console
dir=dirname(rstudioapi::getActiveDocumentContext()$path)     # gets the name of the directory of the active script (works only in R studio)
setwd(dir)

library("R.matlab")
library("reshape2")



################ Specify types of data
prd="funique_control"

################------------------------------------------------ Graph f 

fres=readMat("f_prdData_funique_all.mat")$f[,,1]

#fres=readMat(paste("f_prdData_", prd, ".mat", sep=""))$f[,,1]
fres=data.frame(tank=names(fres), f=matrix(unlist(fres), nrow=length(fres), byrow=T))

fres$tk=colsplit(fres$tank,pattern="\\.", c("tw", "tk", "condition"))[, 2]
fres$condition=colsplit(fres$tank,pattern="\\.", c("tw", "tk", "condition"))[, 3]

fres$condition="BPA0"
fres$study=substr(fres$tk, 1, 5)
fres[fres$tank %in% c("tW.gw124fin", "tW.gw124.BPA100end"), "study"]="gw124end"


library(sciplot)
bargraph.CI(condition, f, data=fres[fres$study=="gw150",], ylim=c(0,1))
bargraph.CI(condition, f, data=fres[fres$study=="gw124",], ylim=c(0,1))
bargraph.CI(condition, f, data=fres[fres$study=="gw124end",], ylim=c(0,1))

temp=fres[which(fres$study %in% c("gw124", "gw150") &
                  fres$condition %in% c("BPA0")),]
bargraph.CI(study, f, data=temp, ylim=c(0,1), leg.lab=c("study 1", "study 2"))


mod=lm(f~condition, fres[fres$study=="gw150",])
anova(mod)
summary(mod)

library(ggplot2)
fres$condition <- factor(fres$condition, levels = c("BPA0", "BPA03", "BPA3", "BPA30", "BPA100", "BPA100end", "BPA300"))

ggplot(data=fres, aes(x=study, y=f, fill=condition)) + 
  geom_bar(position = "dodge", stat = "summary", fun.y = "mean") +                        # Thinner lines
  xlab("studies") + ylab("f") + # Set axis labels+
  # scale_fill_manual(labels=c("control", "BPA0.3", "BPA3", "BPA30", "BPA100", "BPA300"), 
  #                   values=c("#00FF00", "#669900", "#CC3200", "#CB0000", "#650000", "#000000"))+
  scale_fill_manual("legend", values = c("BPA0" = "#00FF00", "BPA03" = "#669900", "BPA3" = "#CC3200", "BPA30"="#CB0000", "BPA100"="#650000", "BPA100end"="#650000", "BPA300"="#000000"))
ggtitle("f for different studies") +     # Set title
  theme_bw()


################------------------------------------------------ Graph diff

# For pred res: 1st column = day, 2nd column = predict, 3rd column = real
predres=readMat(paste("prdData_", prd, ".mat", sep=""))$prdData[,,1]

res_areaError=data.frame(tank=names(predres), areaError=rep(NA, length(names(predres))))

for (tank in names(predres)){
  temp=data.frame(predres[[tank]])
  names(temp)=c("time", "EW", "W")
  temp$diff=(temp$EW-temp$W)/temp$W*100
  eval(parse(text=paste(tank, "=temp", sep="")))
  temp$difft1=c(temp$diff[-1], temp$diff[length(temp$diff)])
  temp$t1=c(temp$time[-1], temp$time[length(temp$time)])
  areaError=sum((temp$t1+temp$time)/2*(temp$diff-temp$difft1))
  res_areaError[res_areaError$tank==tank, "areaError"]=areaError
}

res_areaError$tk=colsplit(res_areaError$tank,pattern="\\.", c("tw", "tk", "condition"))[, 2]
res_areaError$condition=colsplit(res_areaError$tank,pattern="\\.", c("tw", "tk", "condition"))[, 3]

res_areaError[res_areaError$condition=="", "condition"]="BPA0"
res_areaError$study=substr(res_areaError$tk, 1, 5)
res_areaError[res_areaError$tank %in% c("tW.gw124fin", "tW.gw124.BPA100end"), "study"]="gw124end"

res_areaError$condition <- factor(res_areaError$condition, levels = c("BPA0", "BPA03", "BPA3", "BPA30", "BPA100", "BPA100end", "BPA300"))

res_areaError=res_areaError[!res_areaError$study == "gw124end",]
res_areaError$condition <- factor(res_areaError$condition)

ggplot(data=res_areaError, aes(x=study, y=areaError, fill=condition)) + 
  geom_bar(position = "dodge", stat = "summary", fun.y = "mean") +                        # Thinner lines
  xlab("studies") + ylab("AUC for diff") + # Set axis labels+
  # scale_fill_manual(labels=c("control", "BPA0.3", "BPA3", "BPA30", "BPA100", "BPA300"), 
  #                   values=c("#00FF00", "#669900", "#CC3200", "#CB0000", "#650000", "#000000"))+
  scale_fill_manual("legend", values = c("BPA0" = "#00FF00", "BPA03" = "#669900", "BPA3" = "#CC3200", "BPA30"="#CB0000", "BPA100"="#650000", "BPA100end"="#650000", "BPA300"="#000000"))
ggtitle("f for different studies") +     # Set title
theme_bw()




