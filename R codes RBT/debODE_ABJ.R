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
    # ---- calculation of the shape coefficient
    ##-----------------------------------------
    
    if (H < E.Hb){
      s_M = 1
      } else { 
          if (H < E.Hj){                                     
            s_M = L/ Lb
            } else{
              s_M = Lj/ Lb
              }
        }

    ##-----------------------------------------  
    # ---- Correction of p_Am and v by s_M
    ##-----------------------------------------
    
    v = v*s_M             
    p_Am = p_Am*s_M
    
    
    ##-----------------------------------------  
    # ---- Correction of p_Am and v by TC
    ##-----------------------------------------
    
    TC = exp(((T.A)/(T.ref))-((T.A)/(TempC+273.15)))    # Arrhenius correction coeff
    
    v = v*TC               
    p_Am = p_Am*TC
    
    
    ##-----------------------------------------  
    # ---- Growth rate and mobilization rate
    ##-----------------------------------------
    
    # Assimilation 
    if (H>E.Hb){
      pA=f*p_Am*L^2
    } else {
      pA=0
      }
    
    # If p.M varies accross time
    if (length(p.M)>1) {
      pM = approx(c(floor(t),floor(t)+1), 
                  c(p.M[p.M$time==floor(t), 2], p.M[p.M$time==(floor(t)+1), 2]), xout=t)$y
    } else {pM = p.M}
    
    # Mobilization rate (Out of reserves),  J/cm^3
    pC = E/(L^3) * (E.G * v / L + pM) / (kap * E / (L^3) + E.G)      # eq 2.12 book
    
    # Growth rate
    if (kap * pC < pM){ # if p.M higher than energy available
      r = (E * v / (L^4) - pM / kap)    /    (E / (L^3) + E.G * kap_G / kap)     # negativ term because negative numerator
    } else { 
      r = (E * v / (L^4) - pM / kap)    /    (E / (L^3) + E.G / kap)
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
    
    # Extract Lj and Lb dynamically
    if (H <= E.Hb) {
      dLb = max(0, dL)
    } else dLb=0
    
    if (H <= E.Hj) {
      dLj = max(0, dL)
    } else dLj=0
    
    
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