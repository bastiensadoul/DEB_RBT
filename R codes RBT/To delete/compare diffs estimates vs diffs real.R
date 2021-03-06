##################################################################################
##################         Created November 2017          ########################
#####################         by Bastien Sadoul       ############################
##################################################################################
# Calls parameters from .mat for rainbow trout (out of Add my Pet)
# Calls f values for gw124 and gw150, estimated using run_funique_all.m
# Estimates growth using classic deb model and mean f values 
# Estimates growth with varying PARAM and mean f values
# Calculates diff between CLASSIC and VARYING estimates
# Compares with actual diff
##################################################################################

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


#### ------------------------------------
# ---- CALL DEB ODEs
#### ------------------------------------
source(paste(dir,"debODE_ABJ.R", sep="/"))



##########################################################
############## ---- ESTIMATED WEIGHT WITH STANDARD DEB ABJ
##########################################################

dpf=seq(0,1000, by=1)
i=1

for (rep in unique(totreal$variable)) {
  
  study = unique(totreal[totreal$variable == rep, "study"])
  
  #### ------------------------------------
  # ---- Forcing variables
  #### ------------------------------------
  
  f_studies = readMat(gsub("R codes RBT", "/f_prdData_funique_all.mat", dir))$f[,,1]
  
  for_f = names(f_studies)[grep(study, names(f_studies))]
  param_cont$f = mean(unlist(f_studies)[names(f_studies) %in% for_f])

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
  # ---- CALCULATE tj
  #### ------------------------------------
  
  tjcont = LEHovertime_cont[
    which(abs(LEHovertime_cont[,"H"] - param_cont$E.Hj) == min(abs(LEHovertime_cont[,"H"] - param_cont$E.Hj))),"dpf"]
  
  
  #### ------------------------------------
  # ---- SAVE
  #### ------------------------------------
  
  temp_res=data.frame(dpf)
  
  temp_res$estim_W_cont = LEHovertime_cont$estim_W_cont
  mintime = min(totreal[totreal$variable == rep, "dpf"])
  maxtime = max(totreal[totreal$variable == rep, "dpf"])
  
  temp_res$study=study
  temp_res$variable=rep
  
  temp_res$tjcont=tjcont
  
  temp_res = temp_res[which(temp_res$dpf>=mintime & temp_res$dpf<=maxtime),]
  
  if (i==1){
    estim_res_cont = temp_res
  } else {estim_res_cont = rbind(estim_res_cont, temp_res)}
  
  i=i+1
}
  
  

maxtjcont = max(estim_res_cont$tjcont)


##########################################################
############## ---- ESTIMATED WEIGHT WITH VARYING PARAM
##########################################################
### Needs to be in a function to optimize

rm(list=setdiff(ls(), c("totreal", "debODE_ABJ", "param_cont", "estim_res_cont", "dir")))


source(paste(dir, "spring_and_damper_model.R", sep="/"))


LLode <- function(param_spring_damper){

  
# <-
  
  dpf=seq(0,1000, by=1)
  
  #### ------------------------------------
  # ---- PARAMETERS OF THE SPRING AND DAMPER
  #### ------------------------------------
  
  time_for_var=seq(0,length(dpf)+1, by=1)                    # +1 because in debODE_ABJ "floor(t)+1"
  yini = c(0, 0)
  
  # param_spring_damper = data.frame(ks = 2, cs = 200, Fpert_BPA03 = 20, Fpert_BPA3 = 20,
  #                                  Fpert_BPA30 = 30, Fpert_BPA300 = 20, Fpert_BPA100 = 20)
  
  #  # ---- Make sure all spring and damper parameters are above zero
  
  param_spring_damper = abs(param_spring_damper)+0.0001
  
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
  
  for (rep in unique(totreal$variable[totreal$condition != "BPA0"])){
    
    
    print(i)
    
    
    
    param_deb = param_cont
    
    study = unique(totreal[totreal$variable == rep, "study"])
    condition = unique(totreal$condition[totreal$variable == rep])
    eval(parse(text=paste("Fpert <<- Fpert_", condition, sep="")))
    
    
    # ---- Forcing variables
    
    f_studies = readMat(gsub("R codes RBT", "/f_prdData_funique_all.mat", dir))$f[,,1]
    
    for_f = names(f_studies)[grep(study, names(f_studies))]
    param_cont$f = mean(unlist(f_studies)[names(f_studies) %in% for_f])
    
    param_cont$TempC = 8.5   #  en degres C
    
    
    # ---- Initial state
    
    LEH = numeric(6)
    
    LEH[1] = 0.0001     # L
    LEH[2] = 643.562     # E
    LEH[3] = 0   # H
    LEH[4] = 0     # E_R
    LEH[5] = 0     # Lb, we don't know yet. Will be determined by the ode (when L reaches E_Hb)
    LEH[6] = 0     # Lj, we don't know yet. Will be determined by the ode (when L reaches E_Hj)
  
    
    
    # ---- Creates a p_M varying through time
    
    tp_M=ode(y=yini,func=spring_damper_model, times=time_for_var, parms=c(ks,cs), method="ode45")
    tp_M = as.data.frame(tp_M)
    tp_M[, 2] = (tp_M[, 2]+1) * param_deb$p.M
    param_deb$p.M = tp_M[, c(1,2)]
    tvar = tp_M
  
    # ---- ESTIMATED WEIGHT 
    LEHovertime_var = ode(y = LEH, func = debODE_ABJ, times = dpf, 
                           parms = param_deb,
                           method="ode45")
    colnames(LEHovertime_var) = c("dpf", "L", "E", "H", "E_R", "Lb", "Lj")
    
    LEHovertime_var = as.data.frame(LEHovertime_var)
    
    LEHovertime_var$estim_W_var = LEHovertime_var[,"L"]^3 + 
      LEHovertime_var[,"E"] / param_deb$d.E * param_deb$w_E / param_deb$mu.E   # g, wet weight
  
    # ---- Calculate tj
    tjvar = LEHovertime_var[
      which(abs(LEHovertime_var[,"H"] - param_deb$E.Hj) == min(abs(LEHovertime_var[,"H"] - param_deb$E.Hj))),"dpf"]
    
    # ---- SAVE
    temp_res = data.frame(dpf)
    temp_res$estim_W_var = LEHovertime_var$estim_W_var
    mintime = min(totreal[totreal$variable == rep, "dpf"])
    maxtime = max(totreal[totreal$variable == rep, "dpf"])
    
    temp_res$study=study
    temp_res$variable=rep
  
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
  
  # estim_res_cont = data.frame(estim_res_cont)
  # estim_res_var = data.frame(estim_res_var)
  # 
  # estim_res_cont$dpf=as.factor(estim_res_cont$dpf)
  # estim_res_var$dpf=as.factor(estim_res_var$dpf)
  # estim_res_cont$variable=as.factor(estim_res_cont$variable)
  # estim_res_var$variable=as.factor(estim_res_var$variable)
  # estim_res_cont$study=as.factor(estim_res_cont$study)
  # estim_res_var$study=as.factor(estim_res_var$study)
  # estim_res = merge(estim_res_cont, estim_res_var)

  
  estim_res_cont = data.table(estim_res_cont, key = c("dpf", "study", "variable"))
  estim_res_var = data.table(estim_res_var, key = c("dpf", "study", "variable"))
  tot = data.table(totreal,  key = c("dpf", "study", "variable"))
  
  estim_res = merge(estim_res_cont, estim_res_var)

  estim_res$diff_estimates = (estim_res$estim_W_var - estim_res$estim_W_cont) / estim_res$estim_W_cont *100
  
  totfinal = merge(tot, estim_res, all.x=T)

  # Diff before tj
  totfinal = totfinal[totfinal$dpf < maxtjcont,]
  
  diff=totfinal$real_diffW-totfinal$diff_estimates
  diff = diff[!is.na(diff)]
  
  LL =  as.numeric(t(diff) %*% diff)
  return(LL)
}
  
  
#########################################################################################################
####### ---------  OPTIMISATION
#########################################################################################################

# Starting parameters
param_spring_damper = c(ks = 2, cs = 200, Fpert_BPA03 = 1, Fpert_BPA3 = 10,
                                 Fpert_BPA30 = 12, Fpert_BPA300 = 20, Fpert_BPA100 = 15)

# Run the optimization procedure
tmin=0       # start of the pert (when Fpert applies)
tmax=40       # stop of the pert

trace(LLode, exit= quote(print(c(returnValue(), param_spring_damper))))
# totreal = totreal[which(totreal$study == "gw150" & totreal$condition == "BPA30"),]
# param_spring_damper = data.frame(ks = 2, cs = 200,
                                 # Fpert_BPA30 = 12)
MLresults <- optim(param_spring_damper, LLode, method = c("Nelder-Mead"), hessian=F, 
                   control=list(trace=10))
totreal = tot
untrace(LLode)

write.table(unlist(MLresults), paste(dir, "/results_optim/result_optim_", Sys.Date(), ".txt", sep=""), sep = "\t")



########## ALARM WHEN DONE

install.packages("beepr")
library(beepr)
beep()

# Sends txt when script done (works only on mac)
system("osascript -e 'tell application \"Messages\"' -e 'send \"end of script\" to buddy \"moi\"' -e 'end tell'")








#########################################################################################################
####### ---------  PLOT RESULTS
#########################################################################################################
param_spring_damper = MLresults$par

temp = read.table(paste(dir, "/results_optim/result_optim_2017-11-05.txt", sep=""), sep="\t")
param_spring_damper = c(ks = temp[1,1], cs = temp[2,1], Fpert_BPA03 = temp[3,1], Fpert_BPA3 = temp[4,1],
                        Fpert_BPA30 = temp[5,1], Fpert_BPA300 = temp[6,1], Fpert_BPA100 = temp[7,1])

# Run from commented <-

######## ggplot

library(ggplot2)

estim_res$condition = colsplit(estim_res$variable, "_", names=c("a", "condition"))$condition
estim_res=rbind(data.frame(estim_res), data.frame(dpf=dpf, study="gw150", variable="gw150A_BPA0", estim_W_cont=NA, estim_W_var=NA, diff_estimates=0, condition="BPA0"))
estim_res=rbind(data.frame(estim_res), data.frame(dpf=dpf, study="gw124ini", variable="gw124iniA_BPA0", estim_W_cont=NA, estim_W_var=NA, diff_estimates=0, condition="BPA0"))
estim_res=rbind(data.frame(estim_res), data.frame(dpf=dpf, study="gw124fin", variable="gw124fin", estim_W_cont=NA, estim_W_var=NA, diff_estimates=0, condition="BPA0"))



temp=totfinal[totfinal$study=="gw150",]
p = ggplot(temp, aes(x=dpf, y=real_diffW, color=condition)) +
  stat_summary(fun.y = mean, geom = "point", size=5, alpha=0.6) + 
  stat_summary(fun.y = mean, geom = "line", size=1, alpha=0.6) + 
  stat_summary(fun.data = mean_cl_normal, fun.args=list(mult=1), alpha=0.6)+
  geom_line(data=estim_res[which(estim_res$variable %in% c("gw150A_BPA0","gw150A_BPA3","gw150A_BPA30")),], aes(x=dpf, y=diff_estimates, colour=condition), size=2, alpha=1)+
  scale_fill_manual(labels=c("control", "BPA4", "BPA40"), 
                    #    values=c("#00FF00", "#CC3200", "#CB0000"))+
                    values=c("#00FF00", "#CB0000", "#000000"))+
  scale_color_manual(labels=c("control", "BPA4", "BPA40"), 
                     #    values=c("#00FF00", "#CC3200", "#CB0000"))+
                     values=c("#00FF00", "#CB0000", "#000000"))+
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

temp=totfinal[totfinal$study=="gw124ini",]
p = ggplot(temp, aes(x=dpf, y=real_diffW, color=condition)) +
  stat_summary(fun.y = mean, geom = "point", size=5, alpha=0.6) + 
  stat_summary(fun.y = mean, geom = "line", size=1, alpha=0.6) + 
  stat_summary(fun.data = mean_cl_normal, fun.args=list(mult=1), alpha=0.6)+
  geom_line(data=estim_res[which(estim_res$variable %in% c("gw124iniA_BPA0","gw124iniA_BPA03", "gw124ini_BPA3", "gw124ini_BPA30", "gw124iniA_BPA100", "gw124ini_BPA300")),], aes(x=dpf, y=diff_estimates, colour=condition), size=2, alpha=1)+
  scale_fill_manual(labels=c("control", "BPA03", "BPA3", "BPA30", "BPA100", "BPA300"), 
                        values=colorRampPalette(c("green", "red", "black"))(n = 6))+
  scale_color_manual(labels=c("control", "BPA03", "BPA3", "BPA30", "BPA100", "BPA300"), 
                    values=colorRampPalette(c("green", "red", "black"))(n = 6))+
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


temp=totfinal[totfinal$study=="gw124fin",]
p = ggplot(temp, aes(x=dpf, y=real_diffW, color=condition)) +
  stat_summary(fun.y = mean, geom = "point", size=5, alpha=0.6) + 
  stat_summary(fun.y = mean, geom = "line", size=1, alpha=0.6) + 
  stat_summary(fun.data = mean_cl_normal, fun.args=list(mult=1), alpha=0.6)+
  geom_line(data=estim_res[which(estim_res$variable %in% c("gw124fin","gw124finA_BPA03", "gw124finA_BPAend", "gw124fin_BPA3", "gw124fin_BPA30", "gw124fin_BPA300")),], aes(x=dpf, y=diff_estimates, colour=condition), size=2, alpha=1)+
  scale_fill_manual(labels=c("control", "BPA03", "BPA3", "BPA30", "BPA100", "BPA300"), 
                    values=colorRampPalette(c("green", "red", "black"))(n = 6))+
  scale_color_manual(labels=c("control", "BPA03", "BPA3", "BPA30", "BPA100", "BPA300"), 
                     values=colorRampPalette(c("green", "red", "black"))(n = 6))+
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


