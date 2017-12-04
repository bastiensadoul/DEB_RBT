##################################################################################
##################         Created November 2017          ########################
#####################         by Bastien Sadoul       ############################
##################################################################################
# Dif. eq. for a classic deb ABJ model (accel. between birth and metamorphosis)
# BUT: 
#      - f starts at 64dpf. Before f=0
#      - possibility to have f, pM, EG or pAm varying over time
##################################################################################




debODE_ABJ <- function(t, LEH, parms){
  
  with(as.list(parms), {
    
    dLEH=NULL     # Make sure dLEH is null

    
    ##-----------------------------------------
    #------- extract DEB states variables -----
    ##-----------------------------------------
    L       =   LEH[1]
    E     =   LEH[2]
    H       =   LEH[3]
    E.R    =    LEH[4]
    Lb     =   LEH[5]
    Lj     =   LEH[6]
    
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
    
    
    
    ##--------------Options for birth after a given Lb value
    if (acc_after_Lbcont[1] == TRUE & L<acc_after_Lbcont[2]){
      Lb=acc_after_Lbcont[2]
      Lb=as.numeric(Lb)
    }
    if (acc_after_Lbcont[1] == TRUE & L<acc_after_Lbcont[3]){
      Lj=acc_after_Lbcont[3]
      Lj=as.numeric(Lj)
    }
    ##--------------
    
    
    
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
    
    
    
    ##--------------Options for birth after 64dpf or after a given Lb value
    if (acc_after_64dpf == TRUE & t<64){
      s_M=1
    }
    if (acc_after_Lbcont[1] == TRUE & L<acc_after_Lbcont[2]){
      s_M=1
    }
    if (acc_after_Lbcont[1] == TRUE & L>=acc_after_Lbcont[2] & L<acc_after_Lbcont[3]){
      s_M=L/Lb
    }
    ##--------------
    
    
    
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
    
    
    ##--------------Options for birth after 64dpf or after a given Lb value
    if (acc_after_64dpf == TRUE & t<64){
      pA=0
    }
    if (acc_after_Lbcont[1] == TRUE & L<acc_after_Lbcont[2]){
      pA=0
    }
    ##--------------
    
    
    # Mobilization rate (Out of reserves),  J/cm^3
    pC = E/(L^3) * (EG * v / L + pM) / (kap * E / (L^3) + EG)      # eq 2.12 book
    
    # Growth rate
    if (kap * pC < pM){ # if p.M higher than energy available
      r = (E * v / (L^4) - pM / kap)    /    (E / (L^3) + EG * kap_G / kap)     # negativ term because negative numerator
    } else { 
      r = (E * v / (L^4) - pM / kap)    /    (E / (L^3) + EG / kap)
    }
    
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
    
    # Extract Lb dynamically
    if (H <= E.Hb) {
      dLb = max(0, dL)
    } else dLb=0
    
    
    ##--------------Options for birth after 64dpf or after a given Lb value
    if (acc_after_64dpf == TRUE & t<64){
      dLb=max(0, dL)
    }
    if (acc_after_Lbcont[1] == TRUE & L<acc_after_Lbcont[2]){
      dLb=max(0, dL)
    }
    ##--------------
    
    

    # Extract Lj dynamically
    if (H <= E.Hj) {
      dLj = max(0, dL)
    } else dLj=0
    
    
    
    
    ##--------------Options for birth after 64dpf or after a given Lb value
    if (acc_after_64dpf == TRUE & t<64){
      dLj=max(0, dL)
    }
    if (acc_after_Lbcont[1] == TRUE & L<acc_after_Lbcont[3]){
      dLj=max(0, dL)
    }
    ##--------------
    
    
    
    
    ##-----------------------------------------  
    # ---- Return dE, dH and dL
    ##-----------------------------------------

    dLEH[1] <- dL
    dLEH[2] <- dE
    dLEH[3] <- dH
    dLEH[4] <- dE.R
    dLEH[5] <- dLb
    dLEH[6] <- dLj
    list(dLEH)
  })
}
