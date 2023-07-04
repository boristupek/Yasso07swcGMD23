## Corrected SoilR Yasso07Model #######################################################################
#install.packages("SoilR", dependencies = T)
library(SoilR)
#Yasso07Model #original version


#the original SoilR version IS MISSINNG FRACTIONATION TO AWEN!!!!!
#as a result C input is not distributed to A W E N pools but all enters pool A

# FIXED SoilR Yasso07model function is based on Tuomi et al 2011 fortan code and includes:
# 1) AWEN fractionation
# 2) decomopsition dependency on size of litter
# 3) original environmental functions 
# by Boris Tupek boris.tupek@luke.fi
# Feb 2020

#note: 
#1) size WS can be 0, when component litter awens are summed (woody.awen + non.woody.awen)
#2) if function runs on monthly time steps, data input in monthly levels
#   monthly.litter, monthly.mean.temp (deg.C), monthly.precip.sum
#3) daily time stemp is possible if run with daily data, but not tested here, 
#   dailuy.litter*30*12 = upscaled to year, daily.temp, daily.precip*30  = precip upscaled to month)

#Yasso07 model as in Tuomi et al. 2011 ecological applications
Yasso07Modelfi.month <- function (t, #months 1 month 1/12 of year
                            #decomposition rates
                            ksY = c(kA = 0.73, kW = 5.8, kE = 0.29, kN = 0.031, 
                                    kH = 0.0017)/12, #decrease decomposition rates by 12!!!
                            wetlands, # if wetlands = "y" decrease kSY down to 35% (Goll et al. 2015, Kleinen et al. 2021)
                            #transfers and feedbacks
                            pY = c(p1 = 0.48, p2 = 0.01, p3 = 0.83, p4 = 0.99, 
                                   p5 = 0, p6 = 0.01, p7 = 0, p8 = 0, p9 = 0.03, p10 = 0, 
                                   p11 = 0.01, p12 = 0.92, pH = 0.0045),
                            # climate dependence parameters
                            beta1 = 0.096, 
                            beta2 = -0.0014, 
                            gamma = -1.21, 
                            # Woody litter size dependence parameters
                            delta1 = -1.7, 
                            delta2 = 0.86, 
                            r =  -0.306,
                            C0, #initial C , 5 element vector, e.g. c(0,0,0,0,0) 
                            In, # litter C input, data.frame(years, litter), same length as years, if AWEN faractionatoin not provided it has to be in 5 element form for each time step
                            AWEN, #5 element vector, fractionation of plant C litter input to yasso AWEN pools
                            xi = 0, # x1 != 1 will use climate data
                            # xi = 1  will ignore climate data, no climate effect,
                            MT,# MeanTemperature
                            TA, # TemperatureAmplitude = (mothly temp. range)/2
                            PR_mm, # Precipitation_mm
                            WS, # woody size, 0 no effect  for nonwoody, 2 finewoody, 20 coarse woody
                            solver = deSolve.lsoda.wrapper, 
                            pass = FALSE){
  # AWEN fractionation # see https://en.ilmatieteenlaitos.fi/yasso-download-and-support, 
  # Yasso07 user-interface manual (pdf 450 kB) p.13-14
  t_start = min(t)
  t_end = max(t)
  
  
  if(wetlands=="y"){
    ksY = 0.35*ksY
    #return(ksY)
  }
  
  #structural matrix
  Ap = diag(-1, 5, 5)
  Ap[1, 2] = pY[1]
  Ap[1, 3] = pY[2]
  Ap[1, 4] = pY[3]
  Ap[2, 1] = pY[4]
  Ap[2, 3] = pY[5]
  Ap[2, 4] = pY[6]
  Ap[3, 1] = pY[7]
  Ap[3, 2] = pY[8]
  Ap[3, 4] = pY[9]
  Ap[4, 1] = pY[10]
  Ap[4, 2] = pY[11]
  Ap[4, 3] = pY[12]
  Ap[5, 1:4] = pY[13]
  
  # add Woody litter size dependence to structural matrix
  # WS in cm, e.g. 0 for nonwoody, 2 for finewoody, 20 - for coarsewoody litter
  # the effect of wl size as in Eq. 3.1 in model description
  AYS.wsfun <- function(WS){
    ksY.ws = c(ksY[1:4]*(1+delta1*WS+delta2*(WS^2))^(r),ksY[5])
    #print(ksY.ws)
    A1.ws = abs(diag(ksY.ws))
    AYS.ws = Ap %*% A1.ws
    #print(AYS.ws)
    return(AYS.ws)
  }
  AYS.ws <- AYS.wsfun(WS)
  
  #LitterInput
  if (length(In[1,]) == 6){
    LI = as.matrix(In[,2:6]) #first column for years, 2:6 for AWEN
  } else if (length(AWEN) != 5){
    stop("the AWEN litter fractionation 5 element vector must be provided if C litter input is not fractionated to AWEN")
  } else {
    #fractionate liter input to yasso07 A W E N pools
    LA =  matrix(AWEN , nrow=length(t),
                 ncol=5, byrow=TRUE) 
    LI = LA*as.vector(In[,2])  #first column for years
  }
  
  inputFluxes=function(t){
    matrix(nrow = 5, ncol = 1, LI[t,] )
  }
  inputFluxes_tm=BoundInFluxes(
    inputFluxes, 
    t_start,
    t_end)
  
  #environmental effect as function for matrix multiplication
  #this can be replaced by any soil TEMP-SWC function
  y07.Efun.month <- function(MT, PR_mm){
    #MT = Temp
    PR = 12*PR_mm/1000  # *12 convert monthly precipitation to annual sum # convert  from mm to meters 
    
    #monthly temperature and upscaled monthly precipitation to annual level
    TS =  exp(beta1*MT+beta2*MT^2)*(1-exp(gamma*PR)) 
    TS
  }
  xiE=function(t){
    y07.Efun.month(MT[t], PR_mm[t])
  }
  #if xi = 1 replace climate by no effect
  if(xi == 1){
    #woody size but NO environtal effect 
    #if woody size 0 than no woody size too
    AYS.Ews_t=BoundLinDecompOp(
      function(t){AYS.ws},
      t_start,
      t_end
    )
  } else{
    #woody size and environtal effect
    AYS.Ews_t=BoundLinDecompOp(
      function(t){xiE(t)*AYS.ws},
      t_start,
      t_end
    )
  }
  
  #Yasso07 on Carlos Sierra's general model 
  Mod=GeneralModel(t, A=AYS.Ews_t, ivList= C0, inputFluxes =inputFluxes_tm,
                   solver, pass)
  return(Mod)
}

## example #####
yasso.month.example <- function(){
  
##time series
months=seq(from=1,to=10,by=1)#/365)
Litter=data.frame(year=c(1:10),Litter=rnorm(n=10,mean=10,sd=2))
TempData=data.frame(months,Temp=15+sin(2*pi*months)+
                      rnorm(n=length(months),mean=0,sd=1))
j=length(months) # how many months we simulate 
MeanTemperature <- TempData$Temp  
Precipitation <- rep(80,j) # precipitation 80mm #monthly

MT=MeanTemperature
PR_mm=Precipitation
#note conversion from mm to meters in the model's environmental function


# EXAMPLE of fixed model ##
# Modified yasso07 C. Sierra general model WITH environmental effect 
yasso.month <- Yasso07Modelfi.month(months,
                           C0=rep(0,5), #initial carbon
                           AWEN = c(0.52,0.18,0.08,0.2,0), #to separate litter to yasso AWEN pools, this depends on plant organ and species
                           In=Litter,#litter C input (same length as months)
                           xi = 0, # only xi = 1  will replace climate data no climate effect,
                           wetlands = "n", #mineral soils, no further modification of k-rates
                           MT=MT,#MeanTemperature
                           PR_mm=PR_mm,#Precipitation_mm)
                           WS=0) #WS 2 = woody e.g. branches, roots

Ct=getC(yasso.month)
Rt=getReleaseFlux(yasso.month) #respiration

par(mfrow=c(2,1))
#plot carbon pools
matplot( Ct, type="l", lty=1, col=1:5, xlab="months", ylab="Carbon stocks", main ="YM07(example)")
legend("topleft", c("A","W","E","N", "H"), lty=1, col=1:5, bty="n", n = 2)
#plot respiration
matplot(months, Rt, type="l", ylab="Resiration", lty=1, col=1:2)
}

#