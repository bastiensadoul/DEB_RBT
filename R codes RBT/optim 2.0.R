rm(list=ls())
cat("\014")  # To clear the console
dir=dirname(rstudioapi::getActiveDocumentContext()$path)     # gets the name of the directory of the active script (works only in R studio)

library(akima)      # for interpolation in debODE_ABJ.R
library(R.matlab)
library(deSolve)
library(reshape2)

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
##################################################################################
####### ---------   OPTIONS
##################################################################################
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

# Time scale to make estimation on
dpf=seq(0,1069, by=1)

# Spring and damper parameters
param_spring_damper = c(ks = 0.93, cs = 96, Fpert_BPA03 = 0.0001, Fpert_BPA3 = 10,
                        Fpert_BPA30 = 14, Fpert_BPA300 = 28, Fpert_BPA100 = 16)


tmin=0
tmax=40

# Mode of action "p.M", "E.G" or "p_Am"
MoA = "E.G"

# Acceleration? If sM = False --> no acceleration (sM=1)
sM = "TRUE"

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
totreal$study[which(totreal$study=="gw124" & totreal$dpf<=350)]="gw124ini"
totreal$study[which(totreal$study=="gw124" & totreal$dpf>350)]="gw124fin"

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

#--- Calculate f (mean of control)
f_studies = readMat(gsub("R codes RBT", "/f_prdData_funique_all.mat", dir))$f[,,1]

f_gw150 = mean(unlist(f_studies)[names(f_studies) %in% names(f_studies)[grep("gw150", names(f_studies))]])
f_gw124ini = mean(unlist(f_studies)[names(f_studies) %in% names(f_studies)[grep("gw124ini", names(f_studies))]])
f_gw124fin = mean(unlist(f_studies)[names(f_studies) %in% names(f_studies)[grep("gw124fin", names(f_studies))]])

# f linearly interpolated between ini and fin
f_gw124 = data.frame(time = c(1:(length(dpf)+1)))
f_gw124$f=c(rep(f_gw124ini, 350), 
    approx(c(351,438), c(f_gw124ini, f_gw124fin), xout=c(351:438))$y,
    rep(f_gw124fin, length(dpf)+1-438))


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
  
  
  #### ------------------------------------
  # ---- Initial state
  #### ------------------------------------
  
  LEH = numeric(6)
  
  LEH[1] = 0.0001     # L
  LEH[2] = 643.562     # E
  LEH[3] = 0   # H
  LEH[4] = 0     # E_R
  LEH[5] = 0     # Lb, we don't know yet. Will be determined by the ode (when L reaches E_Hb)
  LEH[6] = 0     # Lj, we don't know yet. Will be determined by the ode (when L reaches E_Hj)
  
  
  #### ------------------------------------
  # ---- GETTING ESTIMATED WEIGHTs
  #### ------------------------------------
  
  LEHovertime_cont = ode(y = LEH, func = debODE_ABJ, times = dpf, 
                         parms = param_cont,
                         method="ode23")
  colnames(LEHovertime_cont) = c("dpf", "L", "E", "H", "E_R", "Lb", "Lj")
  
  LEHovertime_cont = as.data.frame(LEHovertime_cont)
  
  LEHovertime_cont$estim_W_cont = LEHovertime_cont[,"L"]^3 + 
    LEHovertime_cont[,"E"] / param_cont$d.E * param_cont$w_E / param_cont$mu.E   # g, wet weight
  
  
  #### ------------------------------------
  # ---- CALCULATE tb and tj
  #### ------------------------------------

  tbcont = LEHovertime_cont[
    which(abs(LEHovertime_cont[,"H"] - param_cont$E.Hb) == min(abs(LEHovertime_cont[,"H"] - param_cont$E.Hb))),"dpf"]
  
  tjcont = LEHovertime_cont[
    which(abs(LEHovertime_cont[,"H"] - param_cont$E.Hj) == min(abs(LEHovertime_cont[,"H"] - param_cont$E.Hj))),"dpf"]
  
  
  #### ------------------------------------
  # ---- SAVE
  #### ------------------------------------
  
  temp_res=data.frame(dpf)
  
  temp_res$estim_W_cont = LEHovertime_cont$estim_W_cont
  temp_res$estim_L_cont = LEHovertime_cont[,"L"]
  mintime = min(totreal[totreal$study2 == repstudy, "dpf"])
  maxtime = max(totreal[totreal$study2 == repstudy, "dpf"])

  temp_res$study2=repstudy

  temp_res$tbcont=tbcont  
  temp_res$tjcont=tjcont
  
  temp_res = temp_res[which(temp_res$dpf>=mintime & temp_res$dpf<=maxtime),]
  
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
                        "param_spring_damper", "tmin", "tmax", "MoA")))


source(paste(dir, "spring_and_damper_model.R", sep="/"))





#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



LLode <- function(param_spring_damper){
  
  
  # <-

  
  #### ------------------------------------
  # ---- PARAMETERS OF THE SPRING AND DAMPER
  #### ------------------------------------
  
  time_for_var=seq(tmin,length(dpf)+1, by=1)                    # +1 because in debODE_ABJ "floor(t)+1"
  yini = c(0, 0)
  
  # param_spring_damper = data.frame(ks = 2, cs = 200, Fpert_BPA03 = 20, Fpert_BPA3 = 20,
  #                                  Fpert_BPA30 = 30, Fpert_BPA300 = 20, Fpert_BPA100 = 20)
  
  #  # ---- Make sure all spring and damper parameters are above zero
  
  param_spring_damper = abs(param_spring_damper)+0.01
  
  #  # ---- Extract parameters
  
  ks <<- param_spring_damper[1]
  # ks = 2      # force of the spring
  cs <<- param_spring_damper[2]
  # cs = 200         # resilience of the damper
  
  Fpert_BPA03 = param_spring_damper[3]
  Fpert_BPA3 = param_spring_damper[4]
  Fpert_BPA30 = param_spring_damper[5]
  Fpert_BPA300 = param_spring_damper[6]
  Fpert_BPA100 = param_spring_damper[7]
  
  #### ------------------------------------
  # ---- LOOP ON ALL TANKS EXCEPT BPA0
  #### ------------------------------------  
  i=1
  
  for (condition_estim in unique(totreal$condition_estim[totreal$condition != "BPA0"])){
    
#    print(i)
    
    
    param_deb = param_cont

    condition = unique(totreal$condition[totreal$condition_estim == condition_estim])
    study2 = unique(totreal$study2[totreal$condition_estim == condition_estim])
    
    eval(parse(text=paste("Fpert <<- Fpert_", condition, sep="")))
    
    
    # ---- Forcing variables
    
    # f specific for each study
    
    if (study2=="gw150"){param_deb$f=f_gw150
    } else {
      param_deb$f=f_gw124
    }
    
    param_deb$TempC = 8.5   #  en degres C
    
    
    # ---- Initial state
    
    LEH = numeric(6)
    
    LEH[1] = 0.0001     # L
    LEH[2] = 643.562     # E
    LEH[3] = 0   # H
    LEH[4] = 0     # E_R
    LEH[5] = 0     # Lb, we don't know yet. Will be determined by the ode (when L reaches E_Hb)
    LEH[6] = 0     # Lj, we don't know yet. Will be determined by the ode (when L reaches E_Hj)
    
    
    
    # ---- Creates a PARAM varying through time
    
    tPARAM=ode(y=yini,func=spring_damper_model, times=time_for_var, parms=c(ks,cs), method="ode45")
    tPARAM = as.data.frame(tPARAM)
    if(MoA=="p_Am"){tPARAM[, 2] = 1/(tPARAM[, 2]+1)-1}
    eval(parse(text= paste("tPARAM[, 2] = (tPARAM[, 2]+1) * param_deb$", MoA, sep="")))
    eval(parse(text= paste("param_deb$", MoA, " = tPARAM[, c(1,2)]", sep="")))

    
    # ---- ESTIMATED WEIGHT 
    LEHovertime_var = ode(y = LEH, func = debODE_ABJ, times = dpf, 
                          parms = param_deb,
                          method="ode45")
    colnames(LEHovertime_var) = c("dpf", "L", "E", "H", "E_R", "Lb", "Lj")
    
    LEHovertime_var = as.data.frame(LEHovertime_var)
    
    LEHovertime_var$estim_W_var = LEHovertime_var[,"L"]^3 + 
      LEHovertime_var[,"E"] / param_deb$d.E * param_deb$w_E / param_deb$mu.E   # g, wet weight
    
    # ---- Calculate tb, tj
    tbvar = LEHovertime_var[
      which(abs(LEHovertime_var[,"H"] - param_deb$E.Hb) == min(abs(LEHovertime_var[,"H"] - param_deb$E.Hb))),"dpf"]
    
    tjvar = LEHovertime_var[
      which(abs(LEHovertime_var[,"H"] - param_deb$E.Hj) == min(abs(LEHovertime_var[,"H"] - param_deb$E.Hj))),"dpf"]
    
    # ---- SAVE
    temp_res = data.frame(dpf)
    temp_res$estim_W_var = LEHovertime_var$estim_W_var
    temp_res$estim_L_var = LEHovertime_var[,"L"]
    mintime = min(totreal[totreal$study2 == study2, "dpf"])
    maxtime = max(totreal[totreal$study2 == study2, "dpf"])
    
    temp_res$study2 = study2
    temp_res$condition_estim = condition_estim
    temp_res$condition = condition
    
    temp_res$tbvar = tbvar
    temp_res$tjvar = tjvar
    
    temp_res = temp_res[which(temp_res$dpf>=mintime & temp_res$dpf<=maxtime),]
    
    if (i==1){
      estim_res_var = temp_res
    } else {estim_res_var = rbind(estim_res_var, temp_res)}
    
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
  
  totfinal = merge(tot, estim_res, 
                   by=c("dpf", "study2", "condition_estim", "condition"),
                   all.x=T)
  
  # Diff before tj
  totfinal = totfinal[totfinal$dpf > totfinal$tjvar,]
  
  diff=totfinal$real_diffW-totfinal$diff_estimates
  diff = diff[!is.na(diff)]
  
  LL =  as.numeric(t(diff) %*% diff)

  
  return(LL)
}


#########################################################################################################
####### ---------  OPTIMISATION
#########################################################################################################

# Starting parameters
param_spring_damper = c(ks = 0.93, cs = 96, Fpert_BPA03 = 0.001, Fpert_BPA3 = 10,
                        Fpert_BPA30 = 14, Fpert_BPA300 = 29, Fpert_BPA100 = 16)

# Run the optimization procedure
tmin=0       # start of the pert (when Fpert applies)
tmax=40       # stop of the pert

trace(LLode, exit= quote(print(c(returnValue(), param_spring_damper))))
# totreal = totreal[which(totreal$study == "gw150" & totreal$condition == "BPA30"),]
# param_spring_damper = data.frame(ks = 2, cs = 200,
# Fpert_BPA30 = 12)
MLresults <- optim(param_spring_damper, LLode, method = c("Nelder-Mead"), hessian=F, 
                   control=list(trace=10, maxit=1000))
totreal = tot
untrace(LLode)

write.table(unlist(MLresults), paste(dir, "/results_optim/result_optim_", Sys.Date(), ".txt", sep=""), sep = "\t")

########## ALARM WHEN DONE


# Sends txt when script done (works only on mac)
system("osascript -e 'tell application \"Messages\"' -e 'send \"end of script\" to buddy \"moi\"' -e 'end tell'")
















#########################################################################################################
####### ---------  PLOTS
#########################################################################################################


# estim_res=rbind(data.frame(estim_res), data.frame(dpf=dpf, study2="gw124", estim_W_var=NA, condition_estim = NA, condition="BPA0", estim_W_cont=NA, tbcont = NA, tjcont = NA, diff_estimates=0))

  
temp=totfinal[which(totfinal$study2=="gw124" & 
                      totfinal$condition %in% c("BPA0", "BPA300")),]
p = ggplot(temp, aes(x=dpf, y=real_diffW, color=condition)) +
  stat_summary(fun.y = mean, geom = "point", size=5, alpha=0.6) + 
  stat_summary(fun.y = mean, geom = "line", size=1, alpha=0.6) + 
  stat_summary(fun.data = mean_cl_normal, fun.args=list(mult=1), alpha=0.6)+
  geom_line(data=estim_res[which(estim_res$condition %in% c("BPA0", "BPA300")),], aes(x=dpf, y=diff_estimates, colour=condition), size=2, alpha=1)+
  scale_fill_manual(labels=c("control", "BPA300"), 
                    values=colorRampPalette(c("green", "black"))(n = 2))+
  scale_color_manual(labels=c("control", "BPA300"), 
                     values=colorRampPalette(c("green", "black"))(n = 2))+
  # geom_point(alpha=0.6,size=5)+
  # geom_line(alpha=0.6,size=1)+
  expand_limits(x=c(0, 1100),y=c(-30,30))+
  #  expand_limits(x=c(0, 600),y=c(-30,30))+
  labs(x="Days post fertilization", y="Mass difference to control (%)") + 
  theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y = element_text(size=16, colour = "black"),
        axis.title.x = element_text(size=16, margin=margin(t=10)), axis.title.y = element_text(size=16, margin=margin(r=10)),
        legend.text = element_text(size=16), legend.title = element_text(size=16),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background=element_rect("grey", fill="white", size=1)
  )
p




eval(parse(text=paste("plot(param_deb$", MoA, "[,1], param_deb$", MoA,"[,2])", sep="")))














