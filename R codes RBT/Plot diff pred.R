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
ODEmethod="euler"
#ODEmethod="ode23"

# Integration step   --> dt<1 might be wrong because approx pb in debODE_ABJ with "floor(t)+1"
dt=1     

# Time scale to make estimation on
dpf=seq(0,1069, by=dt)

# Spring and damper parameters

#param_spring_damper = read.table(paste(dir, "/results_optim/result_optim_p.M_23-nov.-2017 17.27.txt", sep=""), sep = "\t", header=T)
param_spring_damper = read.table(paste(dir, "/results_optim/result_optim_E.G_22-nov.-2017 19.47.txt", sep=""), sep = "\t", header=T)
row.names(param_spring_damper)=substring(row.names(param_spring_damper),5)
param_spring_damper = as.data.frame(t(param_spring_damper))[,c(1:length(t(param_spring_damper)))]
param_spring_damper = unlist(param_spring_damper)

# param_spring_damper["Fpert_BPA300"] = 10

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
E0_gw124 = 643.5622
E0_gw150 = 605.2904


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















##########################################################
############## ---- ESTIMATED WEIGHT WITH VARYING PARAM
##########################################################
### Needs to be in a function to optimize

rm(list=setdiff(ls(), c("totreal", "debODE_ABJ", "param_cont", "estim_res_cont", "dir",
                        "f_gw150", "f_gw124", "dpf",
                        "param_spring_damper", "empirical", "tmin", "tmax", "MoA", "onlyBPA300",
                        "function_var", "studytoplot", "E0_gw124", "E0_gw150", "identical_recovery_time", "ODEmethod",
                        "dt")))


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
# ---- LOOP ON ALL TANKS EXCEPT BPA0
#### ------------------------------------  
i=1

if (onlyBPA300==T){
  toestim = unique(totreal$condition_estim[totreal$condition == "BPA300"])
} else {toestim = unique(totreal$condition_estim[totreal$condition != "BPA0"])}

for (condition_estim in toestim){
  
  # condition_estim = toestim
  print(i)
  
  
  param_deb = param_cont
  
  condition = unique(totreal$condition[totreal$condition_estim == condition_estim])
  study2 = unique(totreal$study2[totreal$condition_estim == condition_estim])
  
  
  # ---- Parameter of the pert model
  eval(parse(text=paste("Fpert <<- Fpert_", condition, sep="")))
  # eval(parse(text=paste("tmax <<- param_spring_damper['tmax_", condition, "']", sep="")))
  if (!is.na(param_spring_damper["tmax"])){tmax <<- param_spring_damper["tmax"]}

  # ---- Forcing variables
  if (study2=="gw150"){param_deb$f=f_gw150
  } else {
    param_deb$f=f_gw124
  }
  
  param_deb$TempC = 8.5   #  en degres C
  
  ### --- Provides E0
  eval(parse(text = paste("E0 = ", "E0_", study2, sep="")))
  
  
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
  stress=ode(y=yini,func=functionforvar, times=time_for_var, parms=c(ks,cs), method=ODEmethod)
  
  stress = as.data.frame(stress)
  tPARAM = stress
  if (function_var == "decreasing_logistic") {tPARAM[, 2] = stress[, 2]*Fpert
  }
  
  # Multiply/Divise stress by coeff
  eval(parse(text=c("iniparam = param_deb$", MoA)))
  if(MoA=="p_Am"){
    stress[, 2] = stress[, 2]/100
    stress[, 2] = -stress[, 2]
  }
  
  tPARAM[, 2] = (stress[, 2]+1) * iniparam
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
  tbvar = LEHovertime_var[
    which(abs(LEHovertime_var[,"H"] - param_deb$E.Hb) == min(abs(LEHovertime_var[,"H"] - param_deb$E.Hb))),"dpf"]
  
  tjvar = LEHovertime_var[
    which(abs(LEHovertime_var[,"H"] - param_deb$E.Hj) == min(abs(LEHovertime_var[,"H"] - param_deb$E.Hj))),"dpf"]
  
  Lbvar = LEHovertime_var[
    which(abs(LEHovertime_var[,"H"] - param_deb$E.Hb) == min(abs(LEHovertime_var[,"H"] - param_deb$E.Hb))),"L"]
  
  Ljvar = LEHovertime_var[
    which(abs(LEHovertime_var[,"H"] - param_deb$E.Hj) == min(abs(LEHovertime_var[,"H"] - param_deb$E.Hj))),"L"]
  
  
  
  # ---- SAVE
  temp_res = data.frame(dpf)
  temp_res$estim_W_var = LEHovertime_var$estim_W_var
  temp_res$estim_L_var = LEHovertime_var[,"L"]
  temp_res$estim_E_var = LEHovertime_var[,"E"]
  temp_res$estim_H_var = LEHovertime_var[,"H"]
  
  mintime = min(totreal[totreal$study2 == study2, "dpf"])
  maxtime = max(totreal[totreal$study2 == study2, "dpf"])
  
  temp_res$study2=study2
  temp_res$condition_estim=condition_estim
  temp_res$condition = condition
  
  temp_res$tbvar = tbvar
  temp_res$tjvar = tjvar
  temp_res$Lbvar = Lbvar
  temp_res$Ljvar = Ljvar
  
  # temp_res = temp_res[which(temp_res$dpf>=mintime & temp_res$dpf<=maxtime),]
  
  names(tPARAM)[2]=MoA
  tPARAM$condition_estim=condition_estim
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

estim_res_cont = data.table(estim_res_cont, key = c("dpf", "study2"))
estim_res_var = data.table(estim_res_var, key = c("dpf", "study2", "condition_estim"))
tot = data.table(totreal,  key = c("dpf", "study2", "condition_estim"))

estim_res = merge(estim_res_var, estim_res_cont)

estim_res$diff_estimates = (estim_res$estim_W_var - estim_res$estim_W_cont) / estim_res$estim_W_cont *100
estim_res$diff_E_estimates = (estim_res$estim_E_var - estim_res$estim_E_cont) / estim_res$estim_E_cont *100
estim_res$diff_L_estimates = (estim_res$estim_L_var - estim_res$estim_L_cont) / estim_res$estim_L_cont *100

estim_res$sMvar = estim_res$Ljvar/estim_res$Lbvar
estim_res$sMcont = estim_res$Ljcont/estim_res$Lbcont

totfinal = merge(tot, estim_res, 
                 by=c("dpf", "study2", "condition_estim", "condition"),
                 all.x=T)

# Calculate sMcont et sMvar
totfinal$sMcont = totfinal$Ljcont/totfinal$Lbcont
totfinal$sMvar = totfinal$Ljvar/totfinal$Lbvar

# # Diff before tj
# totfinal = totfinal[totfinal$dpf > totfinal$tjvar,]

fordiff=totfinal[totfinal$condition!="BPA0",]
fordiff=fordiff[which(!fordiff$dpf %in% c(784.5, 958.5, 616.5)),]
#fordiff=fordiff[fordiff$dpf!=64,]
diff=fordiff$real_diffW-fordiff$diff_estimates
#diff = diff[!is.na(diff)]

LL =  as.numeric(t(diff) %*% diff)


plot(fordiff$dpf, diff)

#########################################################################################################
####### ---------  CALCULATE RE
#########################################################################################################

LL           # Sum of square difference between real diff and estimated diff
# as.numeric(t(diff[-c(1:20)]) %*% diff[-c(1:20)])


#########################################################################################################
####### ---------  CALCULATE CONSEQUENCES
#########################################################################################################

#----------------------  Difference at t=0 of the MoA (equivalent to the max value)

eval(parse(text = paste("estim_param_var$MoA = estim_param_var$", MoA, sep="")))

maxMoA = aggregate(MoA~condition, data=estim_param_var, FUN=max)
maxMoA=rbind(maxMoA, c(condition="BPA0", MoA=iniparam))

maxMoA$MoA=as.numeric(maxMoA$MoA)
maxMoA$concentration = as.numeric(c(0.3, 100, 3, 30, 300, 0))
maxMoA$forcol = factor(maxMoA$concentration, levels = c("0", "0.3", "3", "30", "100", "300"))


#----------------------  Lb value
Lbvar_per_cond = aggregate(Lbvar~condition*study2, data=estim_res_var, FUN=max)

#----------------------  Lj value
Ljvar_per_cond = aggregate(Ljvar~condition*study2, data=estim_res_var, FUN=max)

#----------------------  tb value
tbvar_per_cond = aggregate(tbvar~condition*study2, data=estim_res_var, FUN=unique)

#----------------------  tj value
tjvar_per_cond = aggregate(tjvar~condition*study2, data=estim_res_var, FUN=unique)

#----------------------  merge
sM_table = merge(Lbvar_per_cond, Ljvar_per_cond)
sM_table = merge(sM_table, tbvar_per_cond)
sM_table = merge(sM_table, tjvar_per_cond)
sM_table$sM = sM_table$Ljvar /  sM_table$Lbvar
sM_table = sM_table[order(sM_table$study2, sM_table$condition),]
sM_table = sM_table[c(1,3,4,2,5,6,7),]
t(sM_table)

#########################################################################################################
####### ---------  PLOTS
#########################################################################################################

if (MoA=="E.G"){
  estim_param_var[,2] = estim_param_var[,2]/1000
}


# ------------------------------------------------------------- PLOT Diff pred gw150
  
studytoplot = "gw150"

temp=totfinal[which(totfinal$study2 == studytoplot),]
tempestim = estim_res[estim_res$study2 == studytoplot,]
tempparam=estim_param_var[estim_param_var$condition %in% unique(temp$condition),]
tempparam$var=tempparam[,2]



# Create colour vector
if (studytoplot == "gw150"){
  labvec = c("control", "BPA3", "BPA30")
  colvec=colorRampPalette(c("green", "red", "black"))(n = 6)[c(1,3,4)]
} else {
  labvec = c("control", "BPA0.3", "BPA3", "BPA30", "BPA100", "BPA300")
  colvec=colorRampPalette(c("green", "red", "black"))(n = 6)
  tempestim$condition=factor(tempestim$condition, levels=c("BPA0", "BPA03", "BPA3", "BPA30", "BPA100", "BPA300"))
  temp$condition=factor(temp$condition, levels=c("BPA0", "BPA03", "BPA3", "BPA30", "BPA100", "BPA300"))
  tempparam$condition=factor(tempparam$condition, levels=c("BPA03", "BPA3", "BPA30", "BPA100", "BPA300"))
}

#########  Diff pred
  
#tempestim = tempestim[which(tempestim$dpf>160),]
p_gw150 = ggplot(temp, aes(x=dpf, y=real_diffW, color=condition)) +
  stat_summary(fun.y = mean, geom = "point", size=5, alpha=0.6) + 
  stat_summary(fun.y = mean, geom = "line", size=1, alpha=0.6) + 
  stat_summary(fun.data = mean_cl_normal, fun.args=list(mult=1), alpha=0.6)+
  geom_line(data=tempestim ,
            aes(x=dpf, y=diff_estimates, colour=condition), size=1.5, alpha=1)+
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


#########  Pert

p_gw150_2 = ggplot(data=tempparam, aes(x=time, y=var, col=condition))+
  geom_line(size=1.5)+
  # labs(x="Days post fertilization", y=bquote(.(MoA)*' (KJ.'*~cm^-3*')')) +
  labs(x="Days post fertilization", y=bquote(E[G]*' (KJ.'*~cm^-3*')')) +
  expand_limits(x=c(0, 1100),y=c(0,max(estim_param_var$E.G)))+
  scale_fill_manual(labels=labvec[-1],
                    values=colvec[-1])+
  scale_color_manual(labels=labvec[-1],
                     values=colvec[-1])+
  theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y = element_text(size=16, colour = "black"),
        axis.title.x = element_text(size=16, margin=margin(t=10)), axis.title.y = element_text(size=16, margin=margin(r=10)),
        legend.text = element_text(size=16), legend.title = element_text(size=16),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background=element_rect("grey", fill="white", size=1),
        plot.margin = unit(c(0,3.5,0.5,0.5), "cm"),
        legend.position="none"
  )


# plot_grid(p_gw150, p_gw150_2, nrow=2, rel_heights = c(1/2, 1/3))



# ------------------------------------------------------------- PLOT Diff pred gw124

studytoplot = "gw124"

temp=totfinal[which(totfinal$study2 == studytoplot),]
tempestim = estim_res[estim_res$study2 == studytoplot,]
tempparam=estim_param_var[estim_param_var$condition %in% unique(temp$condition),]
tempparam$var=tempparam[,2]


# Create colour vector
if (studytoplot == "gw150"){
  labvec = c("control", "BPA3", "BPA30")
  colvec=colorRampPalette(c("green", "red", "black"))(n = 6)[c(1,3,4)]
} else {
  labvec = c("control", "BPA0.3", "BPA3", "BPA30", "BPA100", "BPA300")
  colvec=colorRampPalette(c("green", "red", "black"))(n = 6)
  tempestim$condition=factor(tempestim$condition, levels=c("BPA0", "BPA03", "BPA3", "BPA30", "BPA100", "BPA300"))
  temp$condition=factor(temp$condition, levels=c("BPA0", "BPA03", "BPA3", "BPA30", "BPA100", "BPA300"))
  tempparam$condition=factor(tempparam$condition, levels=c("BPA03", "BPA3", "BPA30", "BPA100", "BPA300"))
}

#########  Diff pred

#tempestim = tempestim[which(tempestim$dpf>160),]
p_gw124 = ggplot(temp, aes(x=dpf, y=real_diffW, color=condition)) +
  stat_summary(fun.y = mean, geom = "point", size=5, alpha=0.6) + 
  stat_summary(fun.y = mean, geom = "line", size=1, alpha=0.6) + 
  stat_summary(fun.data = mean_cl_normal, fun.args=list(mult=1), alpha=0.6)+
  geom_line(data=tempestim ,
            aes(x=dpf, y=diff_estimates, colour=condition), size=1.5, alpha=1)+
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


# # Create colour vector
# if (studytoplot == "gw150"){
#   labvec = c("control", "BPA3", "BPA30")
#   colvec=colorRampPalette(c("green", "red", "black"))(n = 6)[c(1,3,4)]
# } else {
#   labvec = c("control", "BPA0.3", "BPA3", "BPA30", "BPA100", "BPA300")
#   colvec=colorRampPalette(c("green", "red", "black"))(n = 6)
#   tempestim$condition=factor(tempestim$condition, levels=c("BPA0", "BPA03", "BPA3", "BPA30", "BPA100", "BPA300"))
#   temp$condition=factor(temp$condition, levels=c("BPA0", "BPA03", "BPA3", "BPA30", "BPA100", "BPA300"))
#   tempparam$condition=factor(tempparam$condition, levels=c("BPA03", "BPA3", "BPA30", "BPA100", "BPA300"))
# }


#########  Pert

p_gw124_2 = ggplot(data=tempparam, aes(x=time, y=var, col=condition))+
  geom_line(size=1.5)+
  labs(x="Days post fertilization", y=bquote(E[G]*' (KJ.'*~cm^-3*')')) +
  expand_limits(x=c(0, 1100),y=c(0,max(estim_param_var$E.G)))+
  scale_fill_manual(labels=labvec[-1],
                    values=colvec[-1])+
  scale_color_manual(labels=labvec[-1],
                     values=colvec[-1])+
  theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y = element_text(size=16, colour = "black"),
        axis.title.x = element_text(size=16, margin=margin(t=10)), axis.title.y = element_text(size=16, margin=margin(r=10)),
        legend.text = element_text(size=16), legend.title = element_text(size=16),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background=element_rect("grey", fill="white", size=1),
        plot.margin = unit(c(0,3.5,0.5,0.5), "cm"),
        legend.position="none"
  )


# plot_grid(p_gw124, p_gw124_2, nrow=2, rel_heights = c(1/2, 1/3))


# ------------------------------------------------------------- PLOT TOTAL Diff pred

p_gw124 <- arrangeGrob(p_gw124, top = textGrob("a", x = unit(0.1, "npc")
                                   , y   = unit(0.8, "npc"), just=c("left","top"),
                                   gp=gpar(col="black", fontsize=30, fontfamily="Times Roman")))
p_gw124_2 <- arrangeGrob(p_gw124_2, top = textGrob("b", x = unit(0.1, "npc")
                                     , y   = unit(1.2, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=30, fontfamily="Times Roman")))
p_gw150 <- arrangeGrob(p_gw150, top = textGrob("c", x = unit(0.1, "npc")
                                     , y   = unit(.8, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=30, fontfamily="Times Roman")))
p_gw150_2 <- arrangeGrob(p_gw150_2, top = textGrob("d", x = unit(0.1, "npc")
                                               , y   = unit(1.2, "npc"), just=c("left","top"),
                                               gp=gpar(col="black", fontsize=30, fontfamily="Times Roman")))




#jpeg(paste(dir, "/Figures/diff_pred.jpg", sep=""), res=600, width=30, height=20, units="cm")
plot_grid(p_gw124,p_gw150, p_gw124_2, p_gw150_2, nrow=2, ncol=2, rel_heights = c(1/2, 1/3, 1/2, 1/3))
#dev.off()



# ------------------------------------------------------------- PLOT Fpert Concentration

maxMoA$MoA = maxMoA$MoA/1000
Fpert = ggplot(data=maxMoA, aes(x=concentration, y=MoA))+
  geom_point(aes(colour=forcol), size=5, alpha=0.6)+ 
  # labs(x="BPA treatment (ng/L)", y=bquote(.(MoA)*' (KJ.'*~cm^-3*')')) +
  labs(x="BPA treatment (ng/L)", y=bquote(E[G]*' (KJ.'*~cm^-3*') at t = 0 dpf')) +
  expand_limits(x=c(0, 300),y=c(0,max(maxMoA$MoA)))+
  geom_smooth(aes(x=concentration, y=MoA),
              method="nls", 
              col="black",
              formula=y~a+Vmax*(1-exp(-x*tau)),
              method.args = list(start=c(a=0,tau=0.05,Vmax=50)),
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
Fpert


#jpeg(paste(dir, "/Figures/EGvsConc.jpg", sep=""), res=600, width=16, height=10, units="cm")
Fpert
#dev.off()




# Additional plots --------------------------------------------------------





######################## CONSEQUENCES ON STATE VARIABLES

#studytoplot = "gw150"
#studytoplot = "gw124"

# temp=estim_res[which(estim_res$study2 == studytoplot),]    # the two studies have diff f
# 
# pW = ggplot(data=temp,
#             aes(x=dpf, y=diff_estimates, colour=condition)) +
#   geom_line(size=2, alpha=1)+
#   expand_limits(x=c(0, 1100),y=c(-30,30))+
#   labs(x="Days post fertilization", y=expression("Mass difference \n to control (%)")) + 
#   theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y = element_text(size=16, colour = "black"),
#         axis.title.x = element_text(size=16, margin=margin(t=10)), axis.title.y = element_text(size=16, margin=margin(r=10)),
#         legend.text = element_text(size=16), legend.title = element_text(size=16),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.background=element_rect("grey", fill="white", size=1),
#         plot.margin = unit(c(0.5,0.5,0.5,1), "cm")
#   )
# 
# 
# pE = ggplot(data=temp,
#             aes(x=dpf, y=diff_E_estimates, colour=condition)) +
#   geom_line(size=2, alpha=1)+
#   expand_limits(x=c(0, 1100),y=c(-30,30))+
#   labs(x="Days post fertilization", y=expression("Reserve difference \n to control (%)")) + 
#   theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y = element_text(size=16, colour = "black"),
#         axis.title.x = element_text(size=16, margin=margin(t=10)), axis.title.y = element_text(size=16, margin=margin(r=10)),
#         legend.text = element_text(size=16), legend.title = element_text(size=16),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.background=element_rect("grey", fill="white", size=1),
#         plot.margin = unit(c(0.5,0.5,0.5,1), "cm")
#   )
# 
# pL = ggplot(data=temp,
#             aes(x=dpf, y=diff_L_estimates, colour=condition)) +
#   geom_line(size=2, alpha=1)+
#   expand_limits(x=c(0, 1100),y=c(-30,30))+
#   labs(x="Days post fertilization", y=expression("Struct. L difference \n to control (%)")) + 
#   theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y = element_text(size=16, colour = "black"),
#         axis.title.x = element_text(size=16, margin=margin(t=10)), axis.title.y = element_text(size=16, margin=margin(r=10)),
#         legend.text = element_text(size=16), legend.title = element_text(size=16),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.background=element_rect("grey", fill="white", size=1),
#         plot.margin = unit(c(0.5,0.5,0.5,1), "cm")
#   )
# 
# plot_grid(pW, pL, pE, nrow=3)




# 
# pWvar = ggplot(data=temp,
#                aes(x=dpf, y=estim_W_var, colour=condition)) +
#   geom_line(size=2, alpha=1)+
#   geom_line(data=unique(temp[,c("dpf", "estim_W_cont")]),
#             aes(x=dpf, y=estim_L_cont), size=1, alpha=1, col="black")+
#   expand_limits(x=c(0, 1100))+
#   labs(x="Days post fertilization", y="Struct Length for PARAM var") + 
#   theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y = element_text(size=16, colour = "black"),
#         axis.title.x = element_text(size=16, margin=margin(t=10)), axis.title.y = element_text(size=16, margin=margin(r=10)),
#         legend.text = element_text(size=16), legend.title = element_text(size=16),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.background=element_rect("grey", fill="white", size=1),
#         plot.margin = unit(c(0.5,0.5,0.5,1), "cm")
#   )
# 
# pLvar = ggplot(data=temp,
#             aes(x=dpf, y=estim_L_var, colour=condition)) +
#   geom_line(size=2, alpha=1)+
#   geom_line(data=unique(temp[,c("dpf", "estim_E_cont")]),
#             aes(x=dpf, y=estim_L_cont), size=1, alpha=1, col="black")+
#   expand_limits(x=c(0, 1100))+
#   labs(x="Days post fertilization", y="Struct Length for PARAM var") + 
#   theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y = element_text(size=16, colour = "black"),
#         axis.title.x = element_text(size=16, margin=margin(t=10)), axis.title.y = element_text(size=16, margin=margin(r=10)),
#         legend.text = element_text(size=16), legend.title = element_text(size=16),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.background=element_rect("grey", fill="white", size=1),
#         plot.margin = unit(c(0.5,0.5,0.5,1), "cm")
#   )
# 
# pEvar = ggplot(data=temp,
#                aes(x=dpf, y=estim_E_var, colour=condition)) +
#   geom_line(size=2, alpha=1)+
#   geom_line(data=unique(temp[,c("dpf", "estim_E_cont")]),
#             aes(x=dpf, y=estim_E_cont), size=1, alpha=1, col="black")+
#   # scale_fill_manual(labels=c("control", "BPA300"), 
#   #                   values=colorRampPalette(c("green", "black"))(n = 2))+
#   # scale_color_manual(labels=c("control", "BPA300"), 
#   #                    values=colorRampPalette(c("green", "black"))(n = 2))+
#   # geom_point(alpha=0.6,size=5)+
#   # geom_line(alpha=0.6,size=1)+
#   expand_limits(x=c(0, 1100))+
#   #  expand_limits(x=c(0, 600),y=c(-30,30))+
#   labs(x="Days post fertilization", y="Struct Length for PARAM var") + 
#   theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y = element_text(size=16, colour = "black"),
#         axis.title.x = element_text(size=16, margin=margin(t=10)), axis.title.y = element_text(size=16, margin=margin(r=10)),
#         legend.text = element_text(size=16), legend.title = element_text(size=16),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.background=element_rect("grey", fill="white", size=1),
#         plot.margin = unit(c(0.5,0.5,0.5,1), "cm")
#   )
# 
# pHvar = ggplot(data=temp,
#                aes(x=dpf, y=estim_H_var, colour=condition)) +
#   geom_line(size=2, alpha=1)+
#   geom_line(data=unique(temp[,c("dpf", "estim_H_cont")]),
#             aes(x=dpf, y=estim_H_cont), size=1, alpha=1, col="black")+
#   # scale_fill_manual(labels=c("control", "BPA300"), 
#   #                   values=colorRampPalette(c("green", "black"))(n = 2))+
#   # scale_color_manual(labels=c("control", "BPA300"), 
#   #                    values=colorRampPalette(c("green", "black"))(n = 2))+
#   # geom_point(alpha=0.6,size=5)+
#   # geom_line(alpha=0.6,size=1)+
#   expand_limits(x=c(0, 1100))+
#   #  expand_limits(x=c(0, 600),y=c(-30,30))+
#   labs(x="Days post fertilization", y="Struct Length for PARAM var") + 
#   theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y = element_text(size=16, colour = "black"),
#         axis.title.x = element_text(size=16, margin=margin(t=10)), axis.title.y = element_text(size=16, margin=margin(r=10)),
#         legend.text = element_text(size=16), legend.title = element_text(size=16),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.background=element_rect("grey", fill="white", size=1),
#         plot.margin = unit(c(0.5,0.5,0.5,1), "cm")
#   )
# 







