##################################################################################
##################         Created November 2017          ########################
#####################         by Bastien Sadoul       ############################
##################################################################################
# Calls parameters from .mat for rainbow trout (out of Add my Pet)
# Calls f values for gw124 and gw150, estimated using run_funique_all.m
# Estimates growth using classic deb model and f values 
# Compares with actual weight values (diff)
##################################################################################

rm(list=ls())
cat("\014")  # To clear the console
dir=dirname(rstudioapi::getActiveDocumentContext()$path)     # gets the name of the directory of the active script (works only in R studio)

library(akima)      # for interpolation in debODE_ABJ.R
library(R.matlab)
library(deSolve)

#########################################################################################################
####### ---------   OPENS REAL DATA
#########################################################################################################

real_tW.gW150    = read.table(paste(dir, "/Data txt/tW.gw150.txt", sep=""), sep='\t', header=T)
real_tW.gW124ini = read.table(paste(dir, "/Data txt/tW.gw124ini.txt", sep=""), sep='\t', header=T)
real_tW.gW124fin = read.table(paste(dir, "/Data txt/tW.gw124fin.txt", sep=""), sep='\t', header=T)

real_tW.gW150$dpf = real_tW.gW150$dpb + 64
real_tW.gW124ini$dpf = real_tW.gW124ini$dpb + 64
real_tW.gW124fin$dpf = real_tW.gW124fin$dpb + 64

real_tW.gW150A = real_tW.gW150[, c(8,2)]
real_tW.gW150B = real_tW.gW150[, c(8,4)]
real_tW.gW150C = real_tW.gW150[, c(8,6)]

real_tW.gW124iniA = real_tW.gW124ini[, c(8,2)]
real_tW.gW124iniB = real_tW.gW124ini[, c(8,4)]
real_tW.gW124iniC = real_tW.gW124ini[, c(8,6)]

real_tW.gW124fin = real_tW.gW124fin[, c(3,2)]

#########################################################################################################
####### ---------   PREPARES THE DEB FUNCTIONS AND ENVIRONMENT
#########################################################################################################


#### ------------------------------------
# ---- CALL DEB PARAMETERS
#### ------------------------------------

#--- Automatically imports from .mat file
param=readMat(gsub("R codes RBT", "Oncorhynchus_mykiss/results_Oncorhynchus_mykiss.mat", dir))$par[,,1]

#--- Add p_Am here to be sure it is the same even when pM changes over time
param$p_Am =  param$z * param$p.M / param$kap  # J/d.cm^2, max assimilation rate (eq coming from the fact that at Lmax, kap*pC=pM. And pC is provided by pA --> kap*pAm=pM)

#--- Missing parameters in .mat
param$w_E = 23.9
param$w_V = 23.9
param$M_V = param$d.V/ param$w_V     # mol/cm^3, volume-specific mass of structure

param$kap_G = param$mu.V * param$M_V / param$E.G     # -, growth efficiency


#### ------------------------------------
# ---- Forcing variables
#### ------------------------------------

f_studies = readMat(gsub("R codes RBT", "/f_prdData_funique_all.mat", dir))$f[,,1]

study = names(f_studies)[1]
param$f = eval(parse(text=paste("f_studies$", study)))

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
####### ---------  COMPARE ESTIMATED GROWTH AND REAL VALUES
#########################################################################################################

dpf=seq(0,2000, by=1)

#### ------------------------------------
# ---- ESTIMATED WEIGHT 
#### ------------------------------------

LEHovertime_cont = ode(y = LEH, func = debODE_ABJ, times = dpf, 
                       parms = param,
                       method="ode23")
colnames(LEHovertime_cont) = c("dpf", "L", "E", "H", "E_R", "Lb", "Lj")

LEHovertime_cont = as.data.frame(LEHovertime_cont)

LEHovertime_cont$estim_W = LEHovertime_cont[,"L"]^3 + 
  LEHovertime_cont[,"E"] / param$d.E * param$w_E / param$mu.E   # g, wet weight

Lphysical = LEHovertime_cont[,"L"]/param$del.M        # cm, physical length

# SAVE
eval(parse(text = paste("estim_", study, "=", "LEHovertime_cont[, c('dpf', 'estim_W')]", sep="")))

#### ------------------------------------
# ---- DIFF REAL weight
#### ------------------------------------

real_tW.gW150A = merge(real_tW.gW150A, estim_tW.gw150A, all.x=T)

real_tW.gW150A$diff_W = (real_tW.gW150A$estim_W - real_tW.gW150A$gw150A) / real_tW.gW150A$gw150A *100
