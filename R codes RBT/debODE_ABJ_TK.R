##################################################################################
##################         Created November 2017          ########################
#####################         by Bastien Sadoul       ############################
##################################################################################
# Dif. eq. for a classic deb ABJ model (accel. between birth and metamorphosis)
# BUT: 
#      - f starts at 64dpf. Before f=0
#      - Stress function as described in 6.28 of Kooijman 2010, but without uptake
#      - Affects MoA (E.G, p.M, or p.Am)
##################################################################################




debODE_ABJ <- function(t, LEHCv, parms){
  
  with(as.list(parms), {
    
    dLEHCv=NULL     # Make sure dLEH is null

    
    # ### Assign for test
    # for (i in c(1:length(parms))){
    #   assign(names(parms)[i], unlist(parms[i]))
    # }
    
    ##-----------------------------------------
    #------- extract DEB states variables -----
    ##-----------------------------------------
    L       =   LEHCv[1]
    E     =   LEHCv[2]
    H       =   LEHCv[3]
    E.R    =    LEHCv[4]
    Lb     =   LEHCv[5]
    Lj     =   LEHCv[6]
    Cv      =    LEHCv[7]
    
    Lm = kap * p_Am / p.M
    
    ##-----------------------------------------
    #------- First feeding occurs at 64dpf -----
    ##-----------------------------------------
    
    if (t<64) {
      f=0
    } else {f=f}
    
    
    ##-----------------------------------------  
    # ---- If f varies over time
    ##-----------------------------------------
    if (length(f)>1) {
      fdt = approx(c(floor(t),floor(t)+1), 
                  c(f[f$time==floor(t), 2], f[f$time==(floor(t)+1), 2]), xout=t)$y
    } else {fdt = f}
    
    
    ##-----------------------------------------  
    # ---- Impact of the stress function on MoA
    ##-----------------------------------------
    # s = (Cv/Ct)                       # 6.42 with no effect conc = 0
    # 
    
    s=0
    
    if (MoA %in% c("E.G", "p_Am", "p.M")){
      eval(parse(text=paste(MoA, " = ", MoA, " * (1+s)", sep="")))  # Correct the param by the stress funct
    }
  
    ##-----------------------------------------  
    # ---- calculation of the shape coefficient
    ##-----------------------------------------
    
    if (H < E.Hb){
      s_M = 1
      } else { 
          if (H < E.Hj){                                     
            s_M = L/ Lb
            #s_M = L/ 0.214616314    # Lb = 0.214616314 for cont
            } else{
              s_M = Lj/ Lb
              #s_M = 0.561039757/ 0.214616314     # Lj = 0.561039757 for cont and f=0.8
              }
        }
    if ("sM" == FALSE){
      s_M=1
    }
    
    if (acc_after_64dpf == TRUE & t<64){
      s_M=1
    }
    
    ##-----------------------------------------  
    # ---- If p.M varies over time
    ##-----------------------------------------
    if (length(p.M)>1) {
      pM = approx(c(floor(t),floor(t)+1), 
                  c(p.M[p.M$time==floor(t), 2], p.M[p.M$time==(floor(t)+1), 2]), xout=t)$y
    } else {pM = p.M}
    
    ##-----------------------------------------  
    # ---- If E.G varies over time
    ##-----------------------------------------
    if (length(E.G)>1) {
      EG = approx(c(floor(t),floor(t)+1), 
                  c(E.G[E.G$time==floor(t), 2], E.G[E.G$time==(floor(t)+1), 2]), xout=t)$y
    } else {EG = E.G}
    
    ##-----------------------------------------  
    # ---- If p.Am varies over time
    ##-----------------------------------------
    if (length(p_Am)>1) {
      pAm = approx(c(floor(t),floor(t)+1), 
                  c(p_Am[p_Am$time==floor(t), 2], p_Am[p_Am$time==(floor(t)+1), 2]), xout=t)$y
    } else {pAm = p_Am}
    
    
    ##-----------------------------------------  
    # ---- Correction of p_Am and v by s_M
    ##-----------------------------------------
    
    v = v*s_M
    pAm = pAm*s_M
    
    
    ##-----------------------------------------  
    # ---- Correction of pM, kJ, p_Am and v by TC
    ##-----------------------------------------
    
    TC = exp(((T.A)/(T.ref))-((T.A)/(TempC+273.15)))    # Arrhenius correction coeff
    
    v = v*TC               
    pAm = pAm*TC
    
    pM = pM*TC               
    k.J = k.J*TC
    
    ##-----------------------------------------  
    # ---- Growth rate and mobilization rate
    ##-----------------------------------------
    
    # Assimilation 
    if (H>E.Hb){
      pA=fdt*pAm*L^2
    } else {
      pA=0
      }
    
    # Mobilization rate (Out of reserves),  J/cm^3
    pC = E/(L^3) * (EG * v / L + pM) / (kap * E / (L^3) + EG)      # eq 2.12 book p. 37
    
    # Growth rate
    if (kap * pC < pM){ # if p.M higher than energy available
      r = (E * v / (L^4) - pM / kap)    /    (E / (L^3) + EG * kap_G / kap)     # negativ term because negative numerator
    } else { 
      r = (E * v / (L^4) - pM / kap)    /    (E / (L^3) + EG / kap)
    }
    
    # print(paste("t = ", t, sep=""))
    # print(paste("kap = ", kap, sep=""))
    # print(paste("Cv = ", Cv, sep=""))
    # print(paste("E.G = ", E.G, sep=""))
    # print(paste("pC = ", pC, sep=""))
    # print(paste("pM = ", pM, sep=""))
    
    ##-----------------------------------------  
    # ---- Generate dE, dH and dL
    ##-----------------------------------------
    
    # Before puberty
    if (H < E.Hp) { 
      dH = (1-kap) * pC * L^3 - k.J * H
    } else { 
      dH = 0
    }
    
    # After puberty
    if (H >= E.Hp) { 
      dE.R = kap.R * ((1-kap) * pC * L^3 - k.J * E.Hp)
    } else { 
      dE.R = 0
    }
        
    # Change in energy in reserve
    dE = pA - pC * (L^3)
    
    # Change in structural length
    dL = L * r / 3 #
    
    # Extract Lj and Lb dynamically
    if (H <= E.Hb) {
      dLb = max(0, dL)
    } else dLb=0
    
    if (H <= E.Hj) {
      dLj = max(0, dL)
    } else dLj=0
    
    
    ##-----------------------------------------  
    # ---- Generate dCv with dilution by growth 6.26 p.223
    ##-----------------------------------------
    
    
    if (t<40) {
      dCv=0
    } else {
      dCv = -Cv * ke * Lm/L - Cv* 3 * dL / L     # Because d(ln(u(x)))/dt = 1/u * du/dt  and 
      # ln(u^3) = 3*ln(u) and
      # l=L/Lm
      
    }

    
    ##-----------------------------------------  
    # ---- Return dE, dH and dL
    ##-----------------------------------------

    dLEHCv[1] <- dL
    dLEHCv[2] <- dE
    dLEHCv[3] <- dH
    dLEHCv[4] <- dE.R
    dLEHCv[5] <- dLb
    dLEHCv[6] <- dLj
    dLEHCv[7] <- dCv
    list(dLEHCv)
  })
}
