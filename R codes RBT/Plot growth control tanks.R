rm(list=ls())
cat("\014")  # To clear the console
dir=dirname(rstudioapi::getActiveDocumentContext()$path)     # gets the name of the directory of the active script (works only in R studio)

library(akima)      # for interpolation in debODE_ABJ.R
library(R.matlab)
library(deSolve)
library(reshape2)
library(ggplot2)
library(data.table)
library(gridExtra)
library(grid)
library(cowplot)

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
##################################################################################
####### ---------   OPTIONS
##################################################################################
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

# ODE method see ?ode (euler is the fastest)
#ODEmethod="euler"
ODEmethod="ode23"

# Integration step   --> dt<1 might be wrong because approx pb in debODE_ABJ with "floor(t)+1"
dt=1     

# Time scale to make estimation on
dpf=seq(0,1069, by=dt)

# Study to plot : "gw150" or "gw124"
studytoplot = "gw124"

# Acceleration? If sM = False --> no acceleration (sM=1)
sM = "TRUE"

# TRUE if acceleration only after f=1 (pM = Lj/Lb after t=64dpf)
acc_after_64dpf = "FALSE"

# Initial reserves
E0_gw124 = 604
E0_gw150 = 644


# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
##################################################################################
####### ---------   REAL VALUES
##################################################################################
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------


#### ------------------------------------
# ---- ALL WEIGHT VALUES IN ONE DATAFRAME
#### ------------------------------------

for (study in c("gw150", "gw124ini", "gw124fin", 
                "gw150_BPA3", "gw150_BPA30",
                "gw124_BPA03", "gw124_BPA3", "gw124_BPA30", 
                "gw124_BPA100", "gw124_BPA100end", "gw124_BPA300")) {
  
  eval(parse(text=paste0(
    "temp = read.table(paste(dir, '/Data txt/tW.", study, ".txt', sep=''), sep='\t', header=T)"
  )))
  
  temp$dpf = temp$dpb + 64
  temp = temp[, -1]
  temp = melt(temp, id = "dpf")
  temp = temp[-grep(temp$variable,pattern = "surv"),]
  names(temp)[names(temp)=="value"] = "gw"
  temp$study = study
  
  if (study == "gw150"){
    totreal = temp
  } else { totreal=rbind(totreal, temp)}
  
}


#### ------------------------------------
# ---- ADD A CONDITION VARIABLE
#### ------------------------------------

totreal$condition = colsplit(totreal$study, "_", names = c("1", "2"))[,2]
totreal$condition[totreal$condition==""]="BPA0"

totreal$study = colsplit(totreal$study, "_", names = c("1", "2"))[,1]
# totreal$study[which(totreal$study=="gw124" & totreal$dpf<=350)]="gw124ini"
# totreal$study[which(totreal$study=="gw124" & totreal$dpf>350)]="gw124fin"
totreal$study[which(totreal$study=="gw124" & totreal$dpf<=391)]="gw124ini"
totreal$study[which(totreal$study=="gw124" & totreal$dpf>391)]="gw124fin"

totreal$variable = as.character(totreal$variable)
totreal$variable[which(totreal$study=="gw124ini")]=
  gsub("gw124", "gw124ini", totreal$variable[which(totreal$study=="gw124ini")])
totreal$variable[which(totreal$study=="gw124fin")]=
  gsub("gw124", "gw124fin", totreal$variable[which(totreal$study=="gw124fin")])
totreal$variable[which(totreal$variable=="gw124finfin")]="gw124fin"
totreal$variable[which(totreal$variable=="gw124iniini")]="gw124ini"

totreal$condition[totreal$condition=="BPA100end"] = "BPA100"

totreal$study2=substring(totreal$study, 1, 5)
totreal$condition_estim = paste(totreal$study2, totreal$condition, sep="_")

totreal$dpf[totreal$dpf==957]=958
totreal$dpf[totreal$dpf==957.5]=958.5

#### ------------------------------------
# ---- CALCULATES A DIFF FOR EACH TANK
#### ------------------------------------

meanBPA0 = aggregate(gw ~ study*dpf, data=totreal[totreal$condition=="BPA0",], FUN=mean)
names(meanBPA0)[3] = "mean_gwBPA0"

totreal = merge(totreal, meanBPA0, all.x=T)

totreal=totreal[which(!is.na(totreal$mean_gwBPA0)),]     # remove those with no mean control
totreal$real_diffW = (totreal$gw - totreal$mean_gwBPA0)/totreal$mean_gwBPA0 * 100


#### ------------------------------------
# ---- GET A RE TANK VERSUS MEAN OF THE TANKS
#### ------------------------------------

mean_diffBPA0 = aggregate(real_diffW ~ study*dpf, data=totreal[totreal$condition=="BPA0",], FUN=mean)
names(mean_diffBPA0)[3] = "mean_diffBPA0"
forRE = merge(totreal[totreal$condition=="BPA0",], mean_diffBPA0)
RE = abs(forRE$mean_diffBPA0 - forRE$real_diffW)/ abs(mean(forRE$real_diffW))
RE = sum(RE)/ length(RE)


# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
##################################################################################
####### ---------   DEB ESTIMATES
##################################################################################
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------


#### ------------------------------------
# ---- CALL DEB PARAMETERS
#### ------------------------------------


#--- Automatically imports from .mat file
param_cont=readMat(gsub("R codes RBT", "Oncorhynchus_mykiss/results_Oncorhynchus_mykiss.mat", dir))$par[,,1]

#--- Add p_Am here to be sure it is the same even when pM changes over time
param_cont$p_Am =  param_cont$z * param_cont$p.M / param_cont$kap  # J/d.cm^2, max assimilation rate (eq coming from the fact that at Lmax, kap*pC=pM. And pC is provided by pA --> kap*pAm=pM)

#--- Missing parameters in .mat
param_cont$w_E = 23.9
param_cont$w_V = 23.9
param_cont$M_V = param_cont$d.V/ param_cont$w_V     # mol/cm^3, volume-specific mass of structure

param_cont$kap_G = param_cont$mu.V * param_cont$M_V / param_cont$E.G     # -, growth efficiency

param_cont$dt=dt

#--- Add options
param_cont$acc_after_64dpf = acc_after_64dpf

#--- Calculate f (mean of control)
f_studies = readMat(gsub("R codes RBT", "/f_prdData_funique_all.mat", dir))$f[,,1]

f_gw150 = mean(unlist(f_studies)[names(f_studies) %in% names(f_studies)[grep("gw150", names(f_studies))]])
f_gw124ini = mean(unlist(f_studies)[names(f_studies) %in% names(f_studies)[grep("gw124ini", names(f_studies))]])
f_gw124fin = mean(unlist(f_studies)[names(f_studies) %in% names(f_studies)[grep("gw124fin", names(f_studies))]])

# f linearly interpolated between ini and fin
f_gw124 = data.frame(time = seq(min(dpf), max(dpf)+1, by = dt))
f_gw124$f = NA
# f_gw124$f[f_gw124$time<=350] = f_gw124ini
# f_gw124$f[f_gw124$time>350 & f_gw124$time<438] = approx(c(350,438), c(f_gw124ini, f_gw124fin), xout=seq(350+dt,438-dt, dt))$y
# f_gw124$f[f_gw124$time>=438] = f_gw124fin
f_gw124$f[f_gw124$time<=391] = f_gw124ini
f_gw124$f[f_gw124$time>391 & f_gw124$time<438] = approx(c(391,438), c(f_gw124ini, f_gw124fin), xout=seq(391+dt,438-dt, dt))$y
f_gw124$f[f_gw124$time>=438] = f_gw124fin

# f_gw150 = 1
# f_gw124 = 1


#### ------------------------------------
# ---- CALL DEB ODEs
#### ------------------------------------
source(paste(dir,"debODE_ABJ.R", sep="/"))


##########################################################
############## ---- ESTIMATED WEIGHT WITH STANDARD DEB ABJ
##########################################################

i=1

for (repstudy in unique(totreal$study2)) {
  
  #### ------------------------------------
  # ---- Forcing variables
  #### ------------------------------------
  
  ### --- f specific to each study
  if (repstudy=="gw150"){param_cont$f=f_gw150
  } else {
    param_cont$f=f_gw124
  }
  
  ### --- Temp
  param_cont$TempC = 8.5   #  en degres C
  
  ### --- Provides E0
  eval(parse(text = paste("E0 = ", "E0_", repstudy, sep="")))
  
  #### ------------------------------------
  # ---- Initial state
  #### ------------------------------------
  
  LEH = numeric(6)
  
  LEH[1] = 0.0001     # L
  LEH[2] = E0     # E
  LEH[3] = 0   # H
  LEH[4] = 0     # E_R
  LEH[5] = 0     # Lb, we don't know yet. Will be determined by the ode (when L reaches E_Hb)
  LEH[6] = 0     # Lj, we don't know yet. Will be determined by the ode (when L reaches E_Hj)
  
  
  #### ------------------------------------
  # ---- GETTING ESTIMATED WEIGHTs
  #### ------------------------------------
  
  LEHovertime_cont = ode(y = LEH, func = debODE_ABJ, times = dpf,
                         parms = param_cont,
                         method=ODEmethod)

  colnames(LEHovertime_cont) = c("dpf", "L", "E", "H", "E_R", "Lb", "Lj")
  
  LEHovertime_cont = as.data.frame(LEHovertime_cont)
  
  LEHovertime_cont$estim_W_cont = LEHovertime_cont[,"L"]^3 + 
    LEHovertime_cont[,"E"] / param_cont$d.E * param_cont$w_E / param_cont$mu.E   # g, wet weight
  
  
  #### ------------------------------------
  # ---- CALCULATE tb and tj, Lb and Lj
  #### ------------------------------------

  tbcont = LEHovertime_cont[
    which(abs(LEHovertime_cont[,"H"] - param_cont$E.Hb) == min(abs(LEHovertime_cont[,"H"] - param_cont$E.Hb))),"dpf"]
  
  tjcont = LEHovertime_cont[
    which(abs(LEHovertime_cont[,"H"] - param_cont$E.Hj) == min(abs(LEHovertime_cont[,"H"] - param_cont$E.Hj))),"dpf"]

  Lbcont = LEHovertime_cont[
    which(abs(LEHovertime_cont[,"H"] - param_cont$E.Hb) == min(abs(LEHovertime_cont[,"H"] - param_cont$E.Hb))),"L"]
    
  Ljcont = LEHovertime_cont[
    which(abs(LEHovertime_cont[,"H"] - param_cont$E.Hj) == min(abs(LEHovertime_cont[,"H"] - param_cont$E.Hj))),"L"]
  
  
  #### ------------------------------------
  # ---- SAVE
  #### ------------------------------------
  
  temp_res=data.frame(dpf)
  
  temp_res$estim_W_cont = LEHovertime_cont$estim_W_cont
  temp_res$estim_L_cont = LEHovertime_cont[,"L"]
  temp_res$estim_E_cont = LEHovertime_cont[,"E"]
  temp_res$estim_H_cont = LEHovertime_cont[,"H"]
  
  mintime = min(totreal[totreal$study2 == repstudy, "dpf"])
  maxtime = max(totreal[totreal$study2 == repstudy, "dpf"])

  temp_res$study2=repstudy

  temp_res$tbcont=tbcont  
  temp_res$tjcont=tjcont
  temp_res$Lbcont=Lbcont
  temp_res$Ljcont=Ljcont
  
  
  # temp_res = temp_res[which(temp_res$dpf>=mintime & temp_res$dpf<=maxtime),]
  # 
  if (i==1){
    estim_res_cont = temp_res
  } else {estim_res_cont = rbind(estim_res_cont, temp_res)}
  
  i=i+1
}


















#########################################################################################################
####### ---------  PLOTS
#########################################################################################################


# estim_res=rbind(data.frame(estim_res), data.frame(dpf=dpf, study2="gw124", estim_W_var=NA, condition_estim = NA, condition="BPA0", estim_W_cont=NA, tbcont = NA, tjcont = NA, diff_estimates=0))

#studytoplot = "gw150"
#studytoplot = "gw124"


studytoplot = "gw124"

temp=totreal[which(totreal$study2 == studytoplot),]
temp = temp[which(temp$condition == "BPA0"),]
tempestim = estim_res_cont[estim_res_cont$study2 == studytoplot,]

p = ggplot(temp, aes(x=dpf, y=gw)) +
  stat_summary(fun.y = mean, geom = "point", size=5, alpha=0.6) + 
  stat_summary(fun.y = mean, geom = "line", size=1, alpha=0.6) + 
  stat_summary(fun.data = mean_cl_normal, fun.args=list(mult=1), alpha=0.6)+
  geom_line(data=tempestim ,
            aes(x=dpf, y=estim_W_cont), size=1, alpha=1, col="red")+
  # scale_fill_manual(labels=c("control", "BPA300"), 
  #                   values=colorRampPalette(c("green", "black"))(n = 2))+
  # scale_color_manual(labels=c("control", "BPA300"), 
  #                    values=colorRampPalette(c("green", "black"))(n = 2))+
  # geom_point(alpha=0.6,size=5)+
  # geom_line(alpha=0.6,size=1)+
  expand_limits(x=c(0, 1100))+
  #  expand_limits(x=c(0, 600),y=c(-30,30))+
  labs(x="Days post fertilization", y="Body mass (g)") + 
  theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y = element_text(size=16, colour = "black"),
        axis.title.x = element_text(size=16, margin=margin(t=10)), axis.title.y = element_text(size=16, margin=margin(r=10)),
        legend.text = element_text(size=16), legend.title = element_text(size=16),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background=element_rect("grey", fill="white", size=1),
        plot.margin = unit(c(0.5,0.5,0.5,1), "cm")
  )

p




studytoplot = "gw150"

temp=totreal[which(totreal$study2 == studytoplot),]
temp = temp[which(temp$condition == "BPA0"),]
tempestim = estim_res_cont[estim_res_cont$study2 == studytoplot,]
tempestim = tempestim[tempestim$dpf<=max(temp$dpf),]

p2 = ggplot(temp, aes(x=dpf, y=gw)) +
  stat_summary(fun.y = mean, geom = "point", size=5, alpha=0.6) + 
  stat_summary(fun.y = mean, geom = "line", size=1, alpha=0.6) + 
  stat_summary(fun.data = mean_cl_normal, fun.args=list(mult=1), alpha=0.6)+
  geom_line(data=tempestim ,
            aes(x=dpf, y=estim_W_cont), size=1, alpha=1, col="red")+
  # scale_fill_manual(labels=c("control", "BPA300"), 
  #                   values=colorRampPalette(c("green", "black"))(n = 2))+
  # scale_color_manual(labels=c("control", "BPA300"), 
  #                    values=colorRampPalette(c("green", "black"))(n = 2))+
  # geom_point(alpha=0.6,size=5)+
  # geom_line(alpha=0.6,size=1)+
  expand_limits(x=c(0, 400))+
  #  expand_limits(x=c(0, 600),y=c(-30,30))+
  labs(x="Days post fertilization", y="Body mass (g)") + 
  theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y = element_text(size=16, colour = "black"),
        axis.title.x = element_text(size=16, margin=margin(t=10)), axis.title.y = element_text(size=16, margin=margin(r=10)),
        legend.text = element_text(size=16), legend.title = element_text(size=16),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background=element_rect("grey", fill="white", size=1),
        plot.margin = unit(c(0.5,0.5,0.5,1), "cm")
  )

p2

p <- arrangeGrob(p, top = textGrob("A", x = unit(0.1, "npc")
                                               , y   = unit(0.9, "npc"), just=c("left","top"),
                                               gp=gpar(col="black", fontsize=30, fontfamily="Times Roman")))
p2 <- arrangeGrob(p2, top = textGrob("B", x = unit(0.1, "npc")
                                   , y   = unit(.9, "npc"), just=c("left","top"),
                                   gp=gpar(col="black", fontsize=30, fontfamily="Times Roman")))


plot_grid(p, p2, ncol=2)















studytoplot = "gw124ini"

temp=totreal[which(totreal$study == studytoplot),]
temp = temp[which(temp$condition == "BPA0"),]
tempestim = estim_res_cont[estim_res_cont$study2 == substr(studytoplot, 1,5),]
tempestim = tempestim[tempestim$dpf<=max(temp$dpf),]

p = ggplot(temp, aes(x=dpf, y=gw)) +
  stat_summary(fun.y = mean, geom = "point", size=5, alpha=0.6) + 
  stat_summary(fun.y = mean, geom = "line", size=1, alpha=0.6) + 
  stat_summary(fun.data = mean_cl_normal, fun.args=list(mult=1), alpha=0.6)+
  geom_line(data=tempestim ,
            aes(x=dpf, y=estim_W_cont), size=1, alpha=1, col="red")+
  # scale_fill_manual(labels=c("control", "BPA300"), 
  #                   values=colorRampPalette(c("green", "black"))(n = 2))+
  # scale_color_manual(labels=c("control", "BPA300"), 
  #                    values=colorRampPalette(c("green", "black"))(n = 2))+
  # geom_point(alpha=0.6,size=5)+
  # geom_line(alpha=0.6,size=1)+
  #expand_limits(x=c(0, 1100))+
  #  expand_limits(x=c(0, 600),y=c(-30,30))+
  labs(x="Days post fertilization", y="Body mass (g)") + 
  theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y = element_text(size=16, colour = "black"),
        axis.title.x = element_text(size=16, margin=margin(t=10)), axis.title.y = element_text(size=16, margin=margin(r=10)),
        legend.text = element_text(size=16), legend.title = element_text(size=16),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background=element_rect("grey", fill="white", size=1),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")
  )

p



studytoplot = "gw124fin"

temp=totreal[which(totreal$study == studytoplot),]
temp = temp[which(temp$condition == "BPA0"),]
tempestim = estim_res_cont[estim_res_cont$study2 == substr(studytoplot, 1,5),]
tempestim = tempestim[tempestim$dpf>=min(temp$dpf),]

p2 = ggplot(temp, aes(x=dpf, y=gw)) +
  stat_summary(fun.y = mean, geom = "point", size=5, alpha=0.6) + 
  stat_summary(fun.y = mean, geom = "line", size=1, alpha=0.6) + 
  stat_summary(fun.data = mean_cl_normal, fun.args=list(mult=1), alpha=0.6)+
  geom_line(data=tempestim ,
            aes(x=dpf, y=estim_W_cont), size=1, alpha=1, col="red")+
  # scale_fill_manual(labels=c("control", "BPA300"), 
  #                   values=colorRampPalette(c("green", "black"))(n = 2))+
  # scale_color_manual(labels=c("control", "BPA300"), 
  #                    values=colorRampPalette(c("green", "black"))(n = 2))+
  # geom_point(alpha=0.6,size=5)+
  # geom_line(alpha=0.6,size=1)+
  #expand_limits(x=c(0, 1100))+
  #  expand_limits(x=c(0, 600),y=c(-30,30))+
  labs(x="Days post fertilization", y="Body mass (g)") + 
  theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y = element_text(size=16, colour = "black"),
        axis.title.x = element_text(size=16, margin=margin(t=10)), axis.title.y = element_text(size=16, margin=margin(r=10)),
        legend.text = element_text(size=16), legend.title = element_text(size=16),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background=element_rect("grey", fill="white", size=1),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")
  )

p2



studytoplot = "gw150"

temp=totreal[which(totreal$study2 == studytoplot),]
temp = temp[which(temp$condition == "BPA0"),]
tempestim = estim_res_cont[estim_res_cont$study2 == studytoplot,]
tempestim = tempestim[tempestim$dpf<=max(temp$dpf),]

p3 = ggplot(temp, aes(x=dpf, y=gw)) +
  stat_summary(fun.y = mean, geom = "point", size=5, alpha=0.6) + 
  stat_summary(fun.y = mean, geom = "line", size=1, alpha=0.6) + 
  stat_summary(fun.data = mean_cl_normal, fun.args=list(mult=1), alpha=0.6)+
  geom_line(data=tempestim ,
            aes(x=dpf, y=estim_W_cont), size=1, alpha=1, col="red")+
  # scale_fill_manual(labels=c("control", "BPA300"), 
  #                   values=colorRampPalette(c("green", "black"))(n = 2))+
  # scale_color_manual(labels=c("control", "BPA300"), 
  #                    values=colorRampPalette(c("green", "black"))(n = 2))+
  # geom_point(alpha=0.6,size=5)+
  # geom_line(alpha=0.6,size=1)+
  expand_limits(x=c(0, 400))+
  #  expand_limits(x=c(0, 600),y=c(-30,30))+
  labs(x="Days post fertilization", y="Body mass (g)") + 
  theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y = element_text(size=16, colour = "black"),
        axis.title.x = element_text(size=16, margin=margin(t=10)), axis.title.y = element_text(size=16, margin=margin(r=10)),
        legend.text = element_text(size=16), legend.title = element_text(size=16),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background=element_rect("grey", fill="white", size=1),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")
  )

p3

p <- arrangeGrob(p, top = textGrob("a", x = unit(0.1, "npc")
                                   , y   = unit(0.8, "npc"), just=c("left","top"),
                                   gp=gpar(col="black", fontsize=30, fontfamily="Times Roman")))
p2 <- arrangeGrob(p2, top = textGrob("b", x = unit(0.1, "npc")
                                     , y   = unit(.8, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=30, fontfamily="Times Roman")))
p3 <- arrangeGrob(p3, top = textGrob("c", x = unit(0.1, "npc")
                                     , y   = unit(.8, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=30, fontfamily="Times Roman")))


jpeg(paste(dir, "/Figures/growth_control.jpg", sep=""), res=600, width=30, height=10, units="cm")
plot_grid(p, p2, p3, ncol=3)
dev.off()








