##################################################################################
##################         Created November 2017          ########################
#####################         by Bastien Sadoul       ############################
##################################################################################
# Calls parameters from .mat for rainbow trout (out of Add my Pet)
# Estimates growth using classic deb model
# Estimates growth using pM, EG or pAm varying over time following a spring and 
# damper model
##################################################################################


rm(list=ls())
cat("\014")  # To clear the console
dir=dirname(rstudioapi::getActiveDocumentContext()$path)     # gets the name of the directory of the active script (works only in R studio)

library(akima)      # for interpolation in debODE_ABJ.R
library(R.matlab)
library(deSolve)


#########################################################################################################
####### ---------   PREPARES THE DEB FUNCTIONS AND ENVIRONMENT
#########################################################################################################


#### ------------------------------------
# ---- CALL DEB PARAMETERS
#### ------------------------------------

#--- Automatically imports from .mat file
param=readMat(gsub("R codes RBT", "Oncorhynchus_mykiss/results_Oncorhynchus_mykiss.mat", dir))$par[,,1]

#--- Calculates p_Am
param$p_Am =  param$z * param$p.M / param$kap  # J/d.cm^2, max assimilation rate (eq coming from the fact that at Lmax, kap*pC=pM. And pC is provided by pA --> kap*pAm=pM)

#--- Missing parameters in .mat
param$w_E = 23.9
param$w_V = 23.9
param$M_V = param$d.V/ param$w_V     # mol/cm^3, volume-specific mass of structure

param$kap_G = param$mu.V * param$M_V / param$E.G     # -, growth efficiency


#### ------------------------------------
# ---- Forcing variables
#### ------------------------------------

param$f = 0.7       # f value for gw124
param$TempC = 8.5   #  en degres C


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
# ---- CALL Differential equations
#### ------------------------------------
source(paste(dir,"debODE_ABJ.R", sep="/"))


#########################################################################################################
####### ---------   COMPARE 'PARAM VARYING' VS 'NOT VARYING'
#########################################################################################################

time=seq(0,2000, by=1)

#### ------------------------------------
# ---- WEIGHT AND LENGTH WHEN PARAM NOT VARYING
#### ------------------------------------

LEHovertime_cont = ode(y = LEH, func = debODE_ABJ, times = time, 
                  parms = param,
                  method="ode23")
colnames(LEHovertime_cont) = c("time", "L", "E", "H", "E_R", "Lb", "Lj")

W = LEHovertime_cont[,"L"]^3 + 
  LEHovertime_cont[,"E"] / param$d.E * param$w_E / param$mu.E   # g, wet weight

Lphysical = LEHovertime_cont[,"L"]/param$del.M        # cm, physical length

# SAVE
Wcont = W
Lcont =  LEHovertime_cont[,"L"]
Lphysicalcont = Lphysical
Econt = LEHovertime_cont[, "E"]
#econt = Econt / (Em * Lcont^3)               
pMcont = param$p.M
Lbcont = LEHovertime_cont[length(LEHovertime_cont[,"Lb"]),"Lb"]
Ljcont = LEHovertime_cont[length(LEHovertime_cont[,"Lj"]),"Lj"]
Hcont=LEHovertime_cont[,"H"]
sMcont = Ljcont/Lbcont
tbcont = LEHovertime_cont[
  which(abs(LEHovertime_cont[,"H"] - param$E.Hb) == min(abs(LEHovertime_cont[,"H"] - param$E.Hb))),"time"]
tjcont = LEHovertime_cont[
  which(abs(LEHovertime_cont[,"H"] - param$E.Hj) == min(abs(LEHovertime_cont[,"H"] - param$E.Hj))),"time"]
tbcont = LEHovertime_cont[
  which(abs(LEHovertime_cont[,"H"] - param$E.Hb) == min(abs(LEHovertime_cont[,"H"] - param$E.Hb))),"time"]

# plot(time, Lcont, main="Structural length over time", type="p", col="green", xlim=c(1,200))
# plot(time, Ljcont, main="Lj over time", type="p", col="green")
# plot(time, Lbcont, main="Lb over time", type="p", col="green")

# xfin=200
# plot(time[c(1:xfin)], Wcont[c(1:xfin)], main="W over time", type="p", col="green")
# plot(time[c(1:xfin)], Lcont[c(1:xfin)]/param$del.M, main="L over time", type="p", col="green")
# plot(time[c(1:xfin)], Econt[c(1:xfin)], main="E over time", type="p", col="green")
# plot(time[c(1:xfin)], Hcont[c(1:xfin)], main="H over time", type="p", col="green")
# abline(v=130.8014)
# abline(h=854.1)
# abline(h=29517*0.99)
# abline(v=1008.6)
# 
# plot(time[c(1:xfin)], econt[c(1:xfin)], main="e over time", type="p", col="green", ylim=c(0,1.1))


# #### ------------------------------------
# # ---- PARAM VARIES FOLLOWING SPRING AND DAMPER
# #### ------------------------------------
#
source(paste(dir, "spring_and_damper_model.R", sep="/"))

time_for_var=seq(0,length(time)+1, by=1)                    # +1 because in debODE_ABJ "floor(t)+1"
yini = c(0, 0)

ks = 2      # force of the spring
cs = 200         # resilience of the damper
Fpert = 20      # force of the perturbation     +- 10 for BPA30,   +- 12 for BPA100    +- 20 for BPA300

tmin=0       # start of the pert (when Fpert applies)
tmax=40       # stop of the pert


#  # ---- Create a p_M varying through time
tp_M=ode(y=yini,func=spring_damper_model, times=time_for_var, parms=c(ks,cs), method="ode45")
tp_M = as.data.frame(tp_M)
tp_M[, 2] = (tp_M[, 2]+1) * param$p.M
param$p.M = tp_M[, c(1,2)]
tvar = tp_M


# ---- Create a E.G varying through time
# tE_G=ode(y=yini,func=spring_damper_model, times=time_for_var, parms=c(ks,cs), method="ode45")
# tE_G = as.data.frame(tE_G)
# tE_G[, 2] = (tE_G[, 2]+1) * param$E.G
# param$E.G = tE_G[, c(1,2)]
# tvar = tE_G

# # ---- Create a p_Am varying through time
# tp_Am=ode(y=yini,func=spring_damper_model, times=time_for_var, parms=c(ks,cs), method="ode45")
# tp_Am = as.data.frame(tp_Am)
# tp_Am[, 2] = (tp_Am[, 2]+1) * param$p_Am
# param$p_Am = tp_Am[, c(1,2)]
# tvar = tp_Am


#### ------------------------------------
# ---- Calculate change on a period of time
#### ------------------------------------

LEHovertime_var = ode(y = LEH, func = debODE_ABJ, times = time,
         parms = param,
         method="ode45")
colnames(LEHovertime_var) = c("time", "L", "E", "H", "E_R", "Lb", "Lj")

W = LEHovertime_var[,"L"]^3 +
  LEHovertime_var[,"E"] / param$d.E * param$w_E / param$mu.E   # g, wet weight

Lphysical = LEHovertime_var[,"L"]/param$del.M        # cm, physical length

# SAVE
Wvar = W
Lvar =  LEHovertime_var[,"L"]
Lphysicalvar = Lphysical
Evar = LEHovertime_var[,"E"]
#evar = Evar / (Em * Lvar^3)
pMvar = param$p.M
Lbvar=LEHovertime_var[length(LEHovertime_var[,"Lb"]),"Lb"]
Ljvar=LEHovertime_var[length(LEHovertime_var[,"Lj"]),"Lj"]
Hvar=LEHovertime_var[,"H"]
sMvar=Ljvar/Lbvar
tbvar = LEHovertime_var[
  which(abs(LEHovertime_var[,"H"] - param$E.Hb) == min(abs(LEHovertime_var[,"H"] - param$E.Hb))),"time"]
tjvar = LEHovertime_var[
  which(abs(LEHovertime_var[,"H"] - param$E.Hj) == min(abs(LEHovertime_var[,"H"] - param$E.Hj))),"time"]


#########################################################################################################
####### ---------   PLOT
#########################################################################################################

diff_W = (Wvar - Wcont)/Wcont*100
diff_pM = (pMvar[pMvar$time %in% time,2]-pMcont)/pMcont *100
diff_E = (Evar-Econt)/Econt *100
diff_L = (Lvar-Lcont)/Lcont * 100


xfin=1000

#plot(time[c(1:xfin)], diff_pM[c(1:xfin)], main="diff p_M varying VS stable p_M over time")
par(mfrow=c(1,1)) 
plot(time_for_var[c(1:xfin)], tvar[c(1:xfin),2], main="PARAM varying over time", type="l", col="red")
abline(v=tjvar)
abline(v=tjcont)
abline(v=tbvar)
abline(v=tbcont)
plot(time[c(1:xfin)], diff_E[c(1:xfin)], main="diff E PARAM varying VS stable PARAM over time", type="p", col="red")
plot(time[c(1:xfin)], diff_L[c(1:xfin)], main="diff L PARAM varying VS stable PARAM over time", type="p", col="red")
plot(time[c(1:xfin)], diff_W[c(1:xfin)], main="diff W PARAM varying VS stable PARAM over time", type="p", col="red")

abline(v=tjvar)
abline(v=tjcont)
abline(v=tbvar)
abline(v=tbcont)
abline(v=64)

# plot(time, evar, main="Scaled reserves over time", type="p", col="red")
# points(time, econt, type="p", col="green")
# 
# plot(time, Evar, main="Reserve over time", type="p", col="red")
# points(time, Econt, type="p", col="green")
# 
plot(time[c(1:xfin)], Wvar[c(1:xfin)], main="Weight over time", type="p", col="red")
points(time[c(1:xfin)], Wcont[c(1:xfin)], type="p", col="green")
# 
# plot(time, Hvar, main="Maturity over time", type="p", col="red")
# points(time, Hcont, type="p", col="green")
# 
# plot(time, LEHovertime_var[,"Lb"], main="Lb over time", type="p", col="red")
# plot(time, LEHovertime_var[,"Lj"], main="Lj over time", type="p", col="red")
# 
xfin=150
plot(time[c(1:xfin)], Lvar[c(1:xfin)], main="Structural length over time", type="l", col="red", ylim=c(0.2, 0.7))
points(time[c(1:xfin)], Lcont[c(1:xfin)], main="Structural length over time", type="l", col="green")
# 
# plot(time, Lphysicalvar, main="Physical length over time", type="p", col="red")
# points(time, Lphysicalcont, main="Physical length over time", type="p", col="green")
# 
# abline(v=tbvar, col="red")
# abline(v=tbcont, col="green")
# 
# abline(v=tjvar, col="red")
# abline(v=tjcont, col="green")

