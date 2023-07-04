
##################################################################################################################
# Reconstruct figures of Tupek et al. 2023: 
# Modeling boreal forest’s mineral soil and peat C stock dynamics with Yasso07 model coupled with updated moisture modifier
#
###################################################################################################################

# Boris Tupek
# Natural Resources Institute Finland (LUKE)
# June 2023
# boris.tupek@luke.fi


## Load calibration functions for SOC and SOC-CO2 ########################################################

# yasso07mioE equilibrium and yasso07mioM monthly time series calibration functions ####################### 
yassomioE.calibr.fun.D12.SOC <- function(calib.params){
  ## Equilibrium SOC simulations with MIO ####
  soc.mio.awenh.e <- data.frame(matrix(NA,9,6))#preallocate the tables
  names(soc.mio.awenh.e) <- c("plot","A", "W", "E","N", "H")
  socmio1.e <- rep(NA,9)
  rh.miolitt1.e <- rep(NA,9)
  for(i in 1:9){
    #i = 1
    #Equilibrium litter
    inawen.wie  <- c(as.numeric(treeunder.w.eq.awen[i,2:5]),0)/10000 #convert from ha/m2 to kg/m2!
    inawen.nwie <- c(as.numeric(treeunder.nw.eq.awen[i,2:5]),0)/10000 #convert from ha/m2 to kg/m2!
    inawen.tie <- inawen.wie+inawen.nwie #total
    In.eq = data.frame(t=1,t(inawen.tie))
    #if(i == 3 |i ==5){ #this doubles the biomass in 3 and 5
    #  inawen.tie <- 1.5*inawen.tie
    #}
    clim.ei <- as.numeric(clim.tn5.swc.e9[i,])
    soiltemp= clim.ei[1]
    soilmoist= clim.ei[2]
    #Equilibrium mio model run
    ymioE <- Yasso07mio.EYM(t.step=0, #equilibrium
                            t=1, #equilibrium
                            C0=rep(0,5), #initial carbon
                            AWEN = rep(0,5), #to separate litter to yasso AWEN pools, this depends on plant organ and species
                            In=In.eq, #litter C input fractionated for AWEN (same length as years)
                            xi = 0, # use climate, only xi = 1  will replace climate data no climate effect,
                            MT= soiltemp,
                            SWC=soilmoist,
                            WS=0, #WS 2 = woody e.g. branches, coarse roots
                            calib.params) 
    
    socE.i <- ymioE[["Ct"]]
    soc.mio.awenh.e[i,] <- c(i,socE.i )
    socmio1.e[i] <- sum(socE.i)
    rh.miolitt1.e[i] <- sum(inawen.tie)
  } 
  ymio.eq.out.names <- c("soc.mio.awenh.e","socmio1.e", "rh.miolitt1.e")
  ymio.eq.out <- list(soc.mio.awenh.e,socmio1.e, rh.miolitt1.e)
  names(ymio.eq.out) <- ymio.eq.out.names 
  
  return(ymio.eq.out)
}

#View(Yasso07mio.EYM)


#calib.params <- ym.mean.params
yassomioM.calibr.fun.D12.SOC.CO2month <- function(calib.params){
  
  #calib.params = refPars.blu$best
  
  ## Equilibrium SOC simulations with MIO ####
  soc.mio.awenh.e_month <- data.frame(matrix(NA,9,6))#preallocate the tables
  names(soc.mio.awenh.e_month) <- c("plot","A", "W", "E","N", "H")
  socmio1.e_month <- rep(NA,9)
  rh.miolitt1.e_month <- rep(NA,9)
  for(i in 1:9){
    #i = 1
    #Equilibrium litter
    inawen.wie  <- c(as.numeric(treeunder.w.eq.awen[i,2:5]),0)/10000 #convert from ha/m2 to kg/m2!
    inawen.nwie <- c(as.numeric(treeunder.nw.eq.awen[i,2:5]),0)/10000 #convert from ha/m2 to kg/m2!
    inawen.tie <- (inawen.wie+inawen.nwie)/12 #total # CONVERT TO MONTH!!! 
    In.eq = data.frame(t=1,t(inawen.tie))
    #if(i == 3 |i ==5){ #this doubles the biomass in 3 and 5
    #  inawen.tie <- 1.5*inawen.tie
    #}
    clim.ei <- as.numeric(clim.tn5.swc.e9[i,])
    soiltemp= clim.ei[1]
    soilmoist= clim.ei[2]
    #Equilibrium mio model run
    ymioE <- Yasso07mio_Month(t.step=0, #equilibrium
                              t=1, #equilibrium
                              C0=rep(0,5), #initial carbon
                              AWEN = rep(0,5), #to separate litter to yasso AWEN pools, this depends on plant organ and species
                              In=In.eq, #litter C input fractionated for AWEN (same length as years)
                              xi = 0, # use climate, only xi = 1  will replace climate data no climate effect,
                              MT= soiltemp,
                              SWC=soilmoist,
                              WS=0, #WS 2 = woody e.g. branches, coarse roots
                              calib.params) 
    
    socE.i <- ymioE[["Ct"]]
    soc.mio.awenh.e_month[i,] <- c(i,socE.i )
    socmio1.e_month[i] <- sum(socE.i)
    rh.miolitt1.e_month[i] <- sum(inawen.tie)
    #}
    ## NOTE. FOR CO2 CALIBRATION TEST!!!
    # strip "_month" extension back to basic name version for soc.mio.awenh.e_month
    # 'soc.mio.awenh.e' = 'soc.mio.awenh.e_month'
    #ymio.eq.out_month.names <- c("soc.mio.awenh.e","socmio1.e_month", "rh.miolitt1.e_month")
    #ymio.eq.out_month <- list(soc.mio.awenh.e_month,socmio1.e_month, rh.miolitt1.e_month)
    #names(ymio.eq.out_month) <- ymio.eq.out_month.names 
    #return(ymio.eq.out_month)
    ##SOC END ########
    
    
    ## Timeseries SOC and CO2 monthly simulations with MIO ####
    #for(i in 1:9){
    #i = 1
    litt.fcm.i <- subset(litter.fc.month.awen_tot, plot== i &  year > 2003) #compatible with climate data 2004:2006
    n.i = 36 #dim(litt.fcm.i)[1]
    months.i = 1:n.i # 3 years 2004:2006 = 36 months #if using days decrease decomposition rates to daily, k/12/30 and use daily litter and temps, swc 
    litt.fcm.i$month.t <- months.i
    litt.fcm.i$H <- 0
    litt.fcm.i[,c("A.tot", "W.tot", "E.tot", "N.tot", "H")] <- litt.fcm.i[,c("A.tot", "W.tot", "E.tot", "N.tot", "H")]/10000 #litter on montly level,convert from kg/ha to kg/m2!
    In.months <- as.matrix(litt.fcm.i[,c("month.t","A.tot", "W.tot", "E.tot", "N.tot", "H")])
    
    #THIS REQUIRES OUTPUT FROM EITHER VERSION (YEARLY, MONTHLY) OF MIO EQUILIBRIUM SOC
    #this requires spinup simulation see yassomioE.calibr.fun.D12.SOC 
    #SAME NAME FROM BOTH soc.mio.awenh.e !!!
    
    #soc.mio.awenh.e <- soc.equilibr.awen 
    
    #output from SOC equilibrium 
    C0.i = as.numeric(soc.mio.awenh.e_month[i,c("A","W","E","N","H")])
    #C0.i = as.numeric(soc.mio.awenh.e[i,c("A","W","E","N","H")]) # initial SOC from moi equilibrium simulations
    
    soiltemp.i = tn5.m[,(i+2)] #+2 to skip year and month columns
    soilmoist.i= sm10.m[,(i+2)]
    
    #Monthly mio model run: note yasso decomposition rates k/12
    ymioM <- Yasso07mio.EYM(t.step=1/12, #equilibrium
                            t=months.i, #equilibrium
                            C0=C0.i, #initial carbon
                            AWEN = rep(0,5), #to separate litter to yasso AWEN pools, this depends on plant organ and species
                            In=In.months, #litter C input fractionated for AWEN (same length as months)
                            xi = 0, # use climate
                            MT= soiltemp.i,
                            SWC=soilmoist.i,
                            WS=0, #WS 2 = woody e.g. branches, coarse roots
                            calib.params) 
    socM.i <- rowSums(ymioM[["Ct"]])
    rhM.i <- rowSums(ymioM[["Rt"]])*44/12 #convert from Kg C m-2 month-1 Kg CO2 m2 month-1
    
    soc.months.i <- cbind(litt.fcm.i[,c("plot", "year", "month")],socM.i)
    rh.months.i <- cbind(litt.fcm.i[,c("plot", "year", "month")],rhM.i)
    
    if(i == 1){
      soc.mio.months.e9 <- soc.months.i
      rh.mio.months.e9 <- rh.months.i
    }else{
      soc.mio.months.e9 <- rbind(soc.mio.months.e9,soc.months.i)
      rh.mio.months.e9 <- rbind(rh.mio.months.e9,rh.months.i)
    }
    
  }
  
  ymio.month.out.names <- c("soc.mio.months.e9","rh.mio.months.e9")
  ymio.month.out <- list(soc.mio.months.e9,rh.mio.months.e9)
  names(ymio.month.out) <- ymio.month.out.names 
  return(ymio.month.out)
}
