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
library(cowplot)
library(grid)

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

# Spring and damper parameters

param_spring_damper = read.table(paste(dir, "/results_optim/result_optim_E.G_22-nov.-2017 15.22.txt", sep=""), sep = "\t", header=T)
row.names(param_spring_damper)=substring(row.names(param_spring_damper),5)
param_spring_damper = as.data.frame(t(param_spring_damper))[,c(1:length(t(param_spring_damper)))]
param_spring_damper = unlist(param_spring_damper)

tmin=0
tmax=0

# Mode of action "p.M", "E.G" or "p_Am"
MoA = "E.G"

# Shall the recovery time be identical (works only )
identical_recovery_time = "TRUE"

# Study to plot : "gw150" or "gw124"
studytoplot = "gw124"


# Acceleration? If sM = False --> no acceleration (sM=1)
sM = "TRUE"

# If only want to test on BPA300
onlyBPA300 = "FALSE"

# TRUE if acceleration only after f=1 (pM = Lj/Lb after t=64dpf)
acc_after_64dpf = "FALSE"

# Choose the function to be used for varying parameter over time 
# ("spring_damper_model", "exp_decrease", "decreasing_logistic", "linearmod")
function_var = "exp_decrease"

# Initial reserves
# E0 = 604   # E0 of gw124
E0 = 644    # E0 of gw150
#E0 = 800 # For tests

# function scaled reserve
# f=0.6954883  # mean of gw124ini
#f=0.5829915  # mean of gw150
f=0.3  # for tests

# Temperature  in C
TempC = 8.5

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


#### ------------------------------------
# ---- CALL DEB ODEs
#### ------------------------------------
source(paste(dir,"debODE_ABJ.R", sep="/"))


##########################################################
############## ---- ESTIMATED WEIGHT WITH STANDARD DEB ABJ
##########################################################

i=1


  #### ------------------------------------
  # ---- Forcing variables
  #### ------------------------------------
  
  ### --- f specific to each study
  param_cont$f=f

  ### --- Temp
  param_cont$TempC = TempC   #  en degres C
  
  ### --- Provides E0
  E0=E0
  
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
  

  estim_res_cont = temp_res




##########################################################
############## ---- ESTIMATED WEIGHT WITH VARYING PARAM
##########################################################
### Needs to be in a function to optimize
# 
# rm(list=setdiff(ls(), c("totreal", "debODE_ABJ", "param_cont", "estim_res_cont", "dir",
#                         "f_gw150", "f_gw124", "dpf",
#                         "param_spring_damper", "empirical", "tmin", "tmax", "MoA", "onlyBPA300",
#                         "function_var", "studytoplot", "E0_gw124", "E0_gw150", "identical_recovery_time", "ODEmethod",
#                         "dt")))


source(paste(dir, "spring_and_damper_model.R", sep="/"))





#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




#### ------------------------------------
# ---- PARAMETERS OF THE SPRING AND DAMPER
#### ------------------------------------


#  # ---- Make sure all spring and damper parameters are above zero

if (function_var == "spring_damper_model"){
  param_spring_damper = abs(param_spring_damper)+0.01
}


#  # ---- Extract parameters

ks <<- param_spring_damper["ks"]
# ks = 2      # force of the spring
cs <<- param_spring_damper["cs"]
# cs = 200         # resilience of the damper

Fpert_BPA03 = abs(param_spring_damper["Fpert_BPA03"])
Fpert_BPA3 = abs(param_spring_damper["Fpert_BPA3"])
Fpert_BPA30 = abs(param_spring_damper["Fpert_BPA30"])
Fpert_BPA300 = abs(param_spring_damper["Fpert_BPA300"])
Fpert_BPA100 = abs(param_spring_damper["Fpert_BPA100"])



#### ------------------------------------
# ---- LOOP ON ALL CONDITIONS EXCEPT BPA0
#### ------------------------------------  
i=1

for (condition in c("BPA03", "BPA3", "BPA30", "BPA300", "BPA100")){
  
  # condition = toestim
  print(i)
  
  
  param_deb = param_cont

  
  
  # ---- Parameter of the pert model
  eval(parse(text=paste("Fpert <<- Fpert_", condition, sep="")))
  # eval(parse(text=paste("tmax <<- param_spring_damper['tmax_", condition, "']", sep="")))
  if (!is.na(param_spring_damper["tmax"])){tmax <<- param_spring_damper["tmax"]}
  
  
  # ---- Forcing variables

  param_deb$f=f
  param_deb$TempC = TempC   #  en degres C
  
  ### --- Provides E0
  E0=E0
  
  
  # ---- Initial state
  
  LEH = numeric(6)
  
  LEH[1] = 0.0001     # L
  LEH[2] = E0     # E
  LEH[3] = 0   # H
  LEH[4] = 0     # E_R
  LEH[5] = 0     # Lb, we don't know yet. Will be determined by the ode (when L reaches E_Hb)
  LEH[6] = 0     # Lj, we don't know yet. Will be determined by the ode (when L reaches E_Hj)
  
  
  
  # ---- Creates a PARAM varying through time
  
  # Initial state of the coefficient multiplying/divising PARAM
  time_for_var=seq(tmin,max(dpf)+dt, by=dt)                    # +dt because in debODE_ABJ "floor(t)+dt"
  if (function_var == "spring_damper_model"){yini = c(0, 0)
  } else if (function_var %in% c("exp_decrease", "linearmod")) {yini = Fpert
  } else if (function_var == "decreasing_logistic") {yini = 0.99
  }
  
  # Selects the function to be used
  eval(parse(text=paste("functionforvar = ", function_var, sep="")))
  
  # Ode for coeff over time
  tPARAM=ode(y=yini,func=functionforvar, times=time_for_var, parms=c(ks,cs), method=ODEmethod)
  
  tPARAM = as.data.frame(tPARAM)
  if (function_var == "decreasing_logistic") {tPARAM[, 2] = tPARAM[, 2]*Fpert
  }
  
  # Multiply/Divise PARAM by coeff
  eval(parse(text=c("iniparam = param_deb$", MoA)))
  if(MoA=="p_Am"){
    tPARAM[, 2] = tPARAM[, 2]/100
    tPARAM[, 2] = -tPARAM[, 2]
  }
  
  tPARAM[, 2] = (tPARAM[, 2]+1) * iniparam
  eval(parse(text= paste("param_deb$", MoA, " = tPARAM[, c(1,2)]", sep="")))
  
  
  # ---- ESTIMATED WEIGHT 
  LEHovertime_var = ode(y = LEH, func = debODE_ABJ, times = dpf, 
                        parms = param_deb,
                        method=ODEmethod)
  colnames(LEHovertime_var) = c("dpf", "L", "E", "H", "E_R", "Lb", "Lj")
  
  LEHovertime_var = as.data.frame(LEHovertime_var)
  
  LEHovertime_var$estim_W_var = LEHovertime_var[,"L"]^3 + 
    LEHovertime_var[,"E"] / param_deb$d.E * param_deb$w_E / param_deb$mu.E   # g, wet weight
  
  
  # ---- Calculate tb, tj, Lb, Lj
  Lbvar = max(LEHovertime_var$Lb)
  Ljvar = max(LEHovertime_var$Lj)
  tbvar = LEHovertime_var[
    which(abs(LEHovertime_var[,"L"] - Lbvar) == min(abs(LEHovertime_var[,"L"] - Lbvar))),"dpf"]
  tjvar = LEHovertime_var[
    which(abs(LEHovertime_var[,"L"] - Ljvar) == min(abs(LEHovertime_var[,"L"] - Ljvar))),"dpf"]
  
  
  
  # ---- SAVE
  temp_res = data.frame(dpf)
  temp_res$estim_W_var = LEHovertime_var$estim_W_var
  temp_res$estim_L_var = LEHovertime_var[,"L"]
  temp_res$estim_E_var = LEHovertime_var[,"E"]
  temp_res$estim_H_var = LEHovertime_var[,"H"]

  temp_res$condition=condition
  temp_res$condition = condition
  
  temp_res$tbvar = tbvar
  temp_res$tjvar = tjvar
  temp_res$Lbvar = Lbvar
  temp_res$Ljvar = Ljvar
  
  # temp_res = temp_res[which(temp_res$dpf>=mintime & temp_res$dpf<=maxtime),]
  
  names(tPARAM)[2]=MoA
  tPARAM$condition=condition
  tPARAM$condition = condition
  
  if (i==1){
    estim_res_var = temp_res
    estim_param_var = tPARAM
  } else {
    estim_res_var = rbind(estim_res_var, temp_res)
    estim_param_var = rbind(estim_param_var, tPARAM)
  }
  
  i=i+1
}

#########################################################################################################
####### ---------  CALCULATES DIFF ESTIMATES AND COMPARE WITH REAL DIFF
#########################################################################################################

### ------ MERGE ALL

estim_res_cont = data.table(estim_res_cont, key = c("dpf"))
estim_res_var = data.table(estim_res_var, key = c("dpf", "condition"))

estim_res = merge(estim_res_var, estim_res_cont)

estim_res$diff_estimates = (estim_res$estim_W_var - estim_res$estim_W_cont) / estim_res$estim_W_cont *100
estim_res$diff_E_estimates = (estim_res$estim_E_var - estim_res$estim_E_cont) / estim_res$estim_E_cont *100
estim_res$diff_L_estimates = (estim_res$estim_L_var - estim_res$estim_L_cont) / estim_res$estim_L_cont *100

estim_res$sMvar = estim_res$Ljvar/estim_res$Lbvar
estim_res$sMcont = estim_res$Ljcont/estim_res$Lbcont


#########################################################################################################
####### ---------  GET MINIMUM VALUE AFTER METAMORPHOSIS OF TREATED
#########################################################################################################


mindiff = aggregate(diff_estimates~condition, estim_res[estim_res$dpf>=estim_res$tjvar,], FUN=min)


#########################################################################################################
####### ---------  PLOTS
#########################################################################################################

#########  Diff pred

labvec = c("BPA0.3", "BPA3", "BPA30", "BPA100", "BPA300")
colvec=colorRampPalette(c("green", "red", "black"))(n = 6)[2:6]

estim_res$condition=factor(estim_res$condition, levels=c("BPA03", "BPA3", "BPA30", "BPA100", "BPA300"))
tempestim=estim_res

p = ggplot(tempestim, aes(x=dpf, y=diff_estimates, color=condition)) +
  geom_line(size=1.5, alpha=1)+
  scale_fill_manual(labels=labvec,
                    values=colvec)+
  scale_color_manual(labels=labvec,
                     values=colvec)+
  # geom_point(alpha=0.6,size=5)+
  # geom_line(alpha=0.6,size=1)+
  expand_limits(x=c(0, 1100),y=c(min(estim_res$diff_estimates),30))+
  #  expand_limits(x=c(0, 600),y=c(-30,30))+
  labs(x="Days post fertilization", y="Mass difference to control (%)") + 
  theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y = element_text(size=16, colour = "black"),
        axis.title.x = element_text(size=16, margin=margin(t=10)), axis.title.y = element_text(size=16, margin=margin(r=10)),
        legend.text = element_text(size=16), legend.title = element_text(size=16),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background=element_rect("grey", fill="white", size=1),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")
  )

p+ggtitle(paste("f = ", f, "     T = ", TempC, "     E0 = ", E0))




#########  Min diff treated

labvec = c("control", "BPA0.3", "BPA3", "BPA30", "BPA100", "BPA300")
colvec=colorRampPalette(c("green", "red", "black"))(n = 6)

mindiff = rbind(mindiff, data.frame(condition="BPA0", diff_estimates="0"))

mindiff$condition=factor(mindiff$condition, levels=c("BPA0", "BPA03", "BPA3", "BPA30", "BPA100", "BPA300"))

mindiff$concentration = as.numeric(c(0.3, 100, 3, 30, 300, 0))
mindiff$diff_estimates = as.numeric(mindiff$diff_estimates)
mindiff$forcol = factor(mindiff$concentration, levels = c("0", "0.3", "3", "30", "100", "300"))


p_mindiff = ggplot(data=mindiff, aes(x=concentration, y=diff_estimates))+
  geom_point(aes(colour=forcol), size=5, alpha=0.6)+ 
  labs(x="BPA treatment (ng/L)", y=bquote("Maximum difference to control (%)")) +
  expand_limits(x=c(0, 300))+
  geom_smooth(aes(x=concentration, y=diff_estimates),
              method="nls", 
              col="black",
              formula=y~-(a+Vmax*(1-exp(-x*tau))),
              method.args = list(start=c(a=5,tau=0.2,Vmax=70)),
              se=FALSE, fullrange=T)+
  scale_fill_manual(labels=labvec,
                    values=colvec)+
  scale_color_manual(labels=labvec,
                     values=colvec)+
  theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y = element_text(size=16, colour = "black"),
        axis.title.x = element_text(size=16, margin=margin(t=10)), axis.title.y = element_text(size=16, margin=margin(r=10)),
        legend.text = element_text(size=16), legend.title = element_text(size=16),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background=element_rect("grey", fill="white", size=1),
        plot.margin = unit(c(0,3.5,0.5,0.5), "cm"),
        legend.position="none"
  )
p_mindiff+ggtitle(paste("f = ", f, "     T = ", TempC, "     E0 = ", E0))


#jpeg(paste(dir, "/Figures/EGvsConc.jpg", sep=""), res=600, width=16, height=10, units="cm")
p_mindiff
#dev.off()
