rm(list=ls())
dir=dirname(rstudioapi::getActiveDocumentContext()$path)     # gets the name of the directory of the current script (works only in R studio)

library(akima)


#########################################################################################################
####### ---------   PREPARE THE DEB FUNCTIONS AND ENVIRONMENT
#########################################################################################################

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
# ---- Forcing variables
#### ------------------------------------

f = 0.8
TempC = 8.5   #  en degres C

#### ------------------------------------
# ---- CALL DEB PARAMETERS
#### ------------------------------------
source(paste(dir,"DEB param RBT.R", sep="/"))

#### ------------------------------------
# ---- CALL Differential equations
#### ------------------------------------
source(paste(dir,"debODE_ABJ.R", sep="/"))


#########################################################################################################
####### ---------   COMPARE 'P_M VARYING' VS 'NOT VARYING'
#########################################################################################################

time=seq(1,1000, by=1)

#### ------------------------------------
# ---- WEIGHT AND LENGTH WHEN P_M NOT VARYING
#### ------------------------------------

LEHovertime_cont = ode(y = LEH, func = debODE_ABJ, times = time, 
                  parms = param,
                  method="ode45")
colnames(LEHovertime_cont) = c("time", "L", "E", "H", "E_R", "Lb", "Lj")

W = LEHovertime_cont[,"L"]^3 + 
  LEHovertime_cont[,"E"] / d_E * w_E / mu_E   # g, wet weight

Lphysical = LEHovertime_cont[,"L"]/del_M        # cm, physical length

# SAVE
Wcont = W
Lcont = Lphysical
pMcont = p_M
Lbcont = LEHovertime_cont[length(LEHovertime_cont[,"Lb"]),"Lb"]
Ljcont = LEHovertime_cont[length(LEHovertime_cont[,"Lj"]),"Lj"]
sMcont = Ljcont/Lbcont

#### ------------------------------------
# ---- P_M VARY FOLLOWING SPRING AND DAMPER
#### ------------------------------------

source(paste(dir, "spring_and_damper_model.R", sep="/"))

time_for_pM=seq(0,length(time)+1, by=1)                    # +1 because in debODE_ABJ "floor(t)+1"
yini = c(0, 0)

ks = 0.00001      # force of the spring
cs = 0.001         # resilience of the damper
Fpert = 300      # force of the perturbation

tmin=0       # start of the pert (when Fpert applies)
tmax=30       # stop of the pert


# Create a p_M varying trought time

tp_M=ode(y=yini,func=spring_damper_model, times=time_for_pM, parms=c(ks,cs), method="ode45")
tp_M = as.data.frame(tp_M)
tp_M[, 2] = tp_M[, 2] + p_M
plot(time_for_pM, tp_M[,2], main="p_M over time", type="l", col="red")
p_M = tp_M[, c(1,2)]

#### ------------------------------------
# ---- Calculate change on a period of time
#### ------------------------------------

LEHovertime_var = ode(y = LEH, func = debODE_ABJ, times = time, 
         parms = param,
         method="ode45")
colnames(LEHovertime_var) = c("time", "L", "E", "H", "E_R", "Lb", "Lj")

W = LEHovertime_var[,"L"]^3 + 
  LEHovertime_var[,"E"] / d_E * w_E / mu_E   # g, wet weight

Lphysical = LEHovertime_var[,"L"]/del_M        # cm, physical length

# SAVE
Wvar = W
Lvar = Lphysical
pMvar = p_M
Lbvar=LEHovertime_var[length(LEHovertime_var[,"Lb"]),"Lb"]
Ljvar=LEHovertime_var[length(LEHovertime_var[,"Lj"]),"Lj"]
sMvar=Ljvar/Lbvar


#########################################################################################################
####### ---------   PLOT
#########################################################################################################

diff_W = (Wvar - Wcont)/Wcont*100
diff_pM = (pMvar[pMvar$time %in% time,2]-pMcont)/pMcont *100


#plot(time, diff_pM, main="diff p_M varying VS stable p_M over time")
plot(time, diff_W, main="diff p_M varying VS stable p_M over time", type="p", col="red")


plot(time, Wvar, main="Weight over time", type="p", col="red")
points(time, Wcont, main="Weight over time", type="p", col="green")



#   
# plot(time, W, main="Wet weight over time")
# plot(time, Lphysical, main="Physical length (cm) over time")
# plot(time, LEHovertime[, 3], main="E over time")
# # plot(time, LEHovertime[, 4], main="H over time")

