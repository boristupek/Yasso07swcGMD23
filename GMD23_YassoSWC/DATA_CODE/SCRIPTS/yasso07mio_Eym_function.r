#########################################################################################################################
## Yasso07mio.EYM model for SOC and CO2 of mineral and organic soils ###################################################################
##
# Boris Tupek boris.tupek@luke.fi
# Natural Resources Institute Finland (LUKE)
# February 2020

# Yasso07mio.EYM_function.r runs in R, R Core Team (2018)(R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/ ).
# The model representation is based on Tuomi et al.(2011) (Tuomi, M., Laiho, R., Repo, A. and Liski, J., 2011. Wood decomposition model for boreal forests. Ecological Modelling, 222(3), pp.709-718),
# and runs on the modeling platform of SoilR package (Sierra et al.2012) including also its own Yasso07Model version (Sierra, C.A., Müller, M. and Trumbore, S.E., 2012. Models of soil organic matter decomposition: the SoilR package, version 1.0. Geoscientific Model Development, 5, pp.1045-1060). 

# Note: BUG : The SoilR Yasso07Model is MISSINNG LITTER FRACTIONATION TO A W E N POOLS!
# The Carbon input in Yasso07Model is not distributed to A W E N pools at a time 1 but all enters pool A and redistributes to pools W E N after. 
# This changes SOC proportions between A W E N pools and model's original decomposition rate.

# Yasso07mio.EYM includes:
# 1) AWEN fractionation
# 2) decomoposition dependency on size of litter
# 3) environmental modifier accounting for reduction of decompositions in water saturated soils 
#    (calibrated for SOC and CO2 of mineral and organic soils Tupek et al 2020 (Biogeosciences))
#    instead of air temperature and precipitation model requires soil temperature and moisture 
#    (potentially useful for higher spatial resolution of soc in the landscape)

# check if SoilR package is installed
if(!require("SoilR",character.only = TRUE)) install.packages("SoilR", dependencies = T)
library("SoilR",character.only = TRUE)

# or for more updated version of the package
# install.packages("devtools", dependencies = T)
# library(devtools)
# devtools::install_github('MPIBGC-TEE/SoilR-exp/pkg', force = T)

## yasso07 model function ################################################################################################
#  for equilibrium SOC (E)
#  and for SOC time series at yearly (Y) or monhtly (M) time steps
##################################################################

# parameters #################
# as in the model description of Tuomi et al. (2011), except for environmental climate parameters
## decomposition rates
 ksY <- data.frame(kA.1 = 0.73, kW.2 = 5.8, kE.3 = 0.29, kN.4 = 0.031, kH.5 = 0.0017)
## transfers and feedbacks
 pY <- data.frame(p1.AW = 0.48, p2.AE = 0.01, p3.AN = 0.83, p4.WA = 0.99, 
                 p5.WE = 0, p6.WN = 0.01, p7.EA = 0, p8.EW = 0, p9.EN = 0.03, p10.NA = 0, 
                 p11.NW = 0.01, p12.NE = 0.92, pH = 0.0045)
## Woody litter size effect parameters
 pWS <- data.frame(delta1 = -1.7, delta2 = 0.86, r = -0.306)
## climate effect parameters
 pE <- data.frame(d=0.999, q=1.743,wopt = 31.69)
 
# vector of all parameters
ym.mean.params <- cbind(ksY,pY,pWS,pE)

# core function #############
Yasso07mio.EYM <- function( t.step, # 0 = Equilibrium SOC, 1 = year or 1/12 = month #t.step months reduces the decomposition rates by 12
                            t, #t= 1 for equilibrium,  or time vector if t.step years or months, note: time vector, litter, soil temperature and moisture data has to be of the same length (n rows)!
                            C0, #initial soil organic carbon , 5 element vector, e.g. c(0,0,0,0,0) 
                            AWEN, #5 element vector, fractionation of plant C litter input to yasso AWEN pools # useful only if running component pools separately otherwise AWEN can be set to c(0,0,0,0,0) 
                            In, # litter C input,  data.frame(col1=timevector, col2 = litter.input.total or cols2:6 = litter.input.fractionated.to.AWEN), same length as time vector (in equilibrium only 1 row), litter either 1 column In[,2] if AWEN is provided, or 5 columns In[,2:6] if AWEN set to c(0,0,0,0,0), if yearly time step litter = yearly sum, if monthly t.step litter = monthly sum (these values should accounts for seasonality in leaves litterfall and finerrot turnover, otherwise it could be annual.litter/12)
                            xi, # xi = 0 uses soil temperature and moisture; xi = 1 sets climate modifier to 1 (ignores climate data) 
                            MT,# Soil Temperature (deg.celsius)
                            SWC, # Soil Water Content (%)
                            WS, # woody size (cm), 0 no effect (can be 0 when litter is given for all forest AWEN components together), otherwise  0 nonwoody (leaves, fineroots),  2 finewoody (branches, coarseroots), 20 coarse woody (snags)
                            parameters ){ #vector of all parameters, same as in the model description of Tuomi et al. 2011,except for environmental climate parameters
  
  # define start and end of the time series
  if(t.step == 0){
    t=1 #for equilibrium only mean litter, temp and moisture is used (no time series)
    }  
    t_start = min(t)
    t_end = max(t)
  
  #yasso07 decomposition rates 
  if(t.step == 0 | t.step == 1){ #for Equilibrium SOC or for Yearly time step original pars from Tuomi et al 2011 
    ksY = as.numeric(parameters[1:5])
    }else if (t.step == 1/12) {# for monthly time step reduce pars by 12 
      ksY = as.numeric(parameters[1:5])/12 } else {
    stop("set t.step to 0 for equilibrium, 1 for years or 1/12 for months")}
  
  #yasso07 transfers and feedbacks
  pY = as.numeric(parameters[6:18])

  #yasso07 wood size litter dependence 
  pWS = as.numeric(parameters[19:21])
  delta1 = pWS[1] 
  delta2 = pWS[2] 
  r = pWS[3] 
  
  #environmental effects of soil temperature and moisture 
  pE = as.numeric(parameters[22:24])
  d= pE[1]
  q= pE[2]
  wopt = pE[3]
  
  #structural matrix as in the model description of Tuomi et al. 2011 
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
  
  # structural matrix with decomposition rates and their dependence on the woody litter size
  # the effect of wl size as in Eq. 3.1 in model description of Tuomi et al. 2011 
  AYS.wsfun <- function(WS){
    ksY.ws = c(ksY[1:4]*(1+delta1*WS+delta2*(WS^2))^(r),ksY[5])
    A1.ws = abs(diag(ksY.ws))
    AYS.ws = Ap %*% A1.ws
    
    return(AYS.ws)
  }
  AYS.ws <- AYS.wsfun(WS)
  
  #environmental effect as function for matrix multiplication
  # based on Davidson et al. (2012).
  y07.EMfun <- function(MT,SWC){ 
    EM= d^((wopt - SWC)^2)*q^((MT-10)/10)
    EM #environmental modifier
  }
  xiE=function(t){
    y07.EMfun(MT[t],SWC[t])
  }
  #if xi = 1 replace climate modifier by 1 no effect
  if(xi == 1){
    #NO environtal effect 
    #if woody size 0 than no woody size is also ignored
    #
    #structural matrix for equilibrium
    AYS.wsE <- 1*AYS.ws
    #matrix for time series
    AYS.Ews_t=BoundLinDecompOp(
      function(t){AYS.ws},
      t_start,
      t_end
    )
  } else{
    #woody size and environtal effect
    #
    #structural matrix for equilibrium
    AYS.wsE <- xiE(1)*AYS.ws
    #matrix for time series
    AYS.Ews_t=BoundLinDecompOp(
      function(t){xiE(t)*AYS.ws},
      t_start,
      t_end
    )
  }
  
  #LitterInput
  if (length(In[1,]) == 6){
    LI = as.matrix(In[,2:6]) #first column for years, 2:6 for AWEN
  } else if (length(AWEN) != 5){
    stop("the AWEN fractionation  must be provided as 5 element vector, AWEN can be set to zero c(0,0,0,0,0) if C litter input is fractionated; refer to AWEN fractions of litter in Yasso07 user-interface manual (pdf 450kB) p.13-14
           \nsee https://en.ilmatieteenlaitos.fi/yasso-download-and-support")
  } else {
    #fractionate liter input to yasso07 A W E N pools
    LA =  matrix(AWEN , nrow=length(t),
                 ncol=5, byrow=TRUE) 
    LI = LA*as.vector(In[,2])  #first column for years
  }
  
  # return SOC in equilibrium 
  if(t.step==0){
        #equilibrium soc solve analyticaly Sierra et al. (2018) : Sierra, C.A., Ceballos-Núñez, V., Metzler, H. and Müller, M., 2018. Representing and understanding the carbon cycle using the theory of compartmental dynamical systems. Journal of Advances in Modeling Earth Systems, 10(8), pp.1729-1734.
        #xss = -B^{-1}(t)*u(t)
        socE = -1*solve(AYS.wsE)%*%as.numeric(LI)
        socE = as.numeric(socE)
        
        Esoc.rh <- list(socE, as.numeric(LI))
        names(Esoc.rh) <-c("Ct", "Rt")
    return(Esoc.rh)  #equilibrium SOC of A W E N H pools and  respiration at equilibrium = litter
    
  } else{ 
  #or return Model object of SOC and CO2 timeseries
  inputFluxes=function(t){
    matrix(nrow = 5, ncol = 1, LI[t,] )
    }
  inputFluxes_tm=BoundInFluxes(
    inputFluxes, 
    t_start,
    t_end)
  
  #Yasso07mio.EMY runs on platform of SoilR general model Sierra et al. (2012) 
  Mod=GeneralModel(t, A=AYS.Ews_t, ivList= C0, inputFluxes =inputFluxes_tm,
                   solver=deSolve.lsoda.wrapper, pass = FALSE)
  
  Ct=data.frame(getC(Mod))
  Rt=data.frame(getReleaseFlux(Mod)) #respiration
  
  names(Ct)<-c("A","W","E","N","H")
  names(Rt)<-c("A","W","E","N","H")
  
  CRt <-list(Ct,Rt)
  names(CRt)<-c("Ct","Rt")
  
  return(CRt) #timeseries SOC and CO2 in each  A W E N H pool
  }
  
}



############################################################################################################
## FEW EXAMPLES ############################################################################################

Yasso07mio.EYM.examples <- function(){
  
#read the parameters and the model function from above
#ym.mean.params <- cbind(ksY,pY,pWS,pE)
  
## eqiulibrium
#annual mean
SoilTemperature <- 6 #deg.C
SWC.mineral<- 15 #%
SWC.organic<- 60 #%
Litter.eq = data.frame(time=1, litter = 5 ) #t C ha 
climate.eq = data.frame(SoilTemperature=6,SWC.mineral=15,SWC.organic=60)

print("Equilibrium data:")
print(c(climate.eq,Litter.eq))
AWEN = c(0.52,0.18,0.08,0.2,0)
print(c("AWEN=",AWEN))

for(i in 1:2){
 # i = 1
SWC <- c(SWC.mineral,SWC.organic)[i]
  
yassomio.eq <- Yasso07mio.EYM(t.step=0, #eqiulibrium
                                t=1, #eqiulibrium
                             C0=rep(0,5), #initial carbon
                             AWEN = c(0.52,0.18,0.08,0.2,0), #to separate litter to yasso AWEN pools, this depends on plant organ and species
                             In=Litter.eq,#litter C input (same length as years)
                             xi = 0, # only xi = 1  will replace climate data no climate effect,
                             MT= SoilTemperature,
                             SWC=SWC,
                             WS=0, #WS 2 = woody e.g. branches, coarse roots
                             ym.mean.params) 
socE <- yassomio.eq[["Ct"]]
#names(socE) <- c("A","W","E","N","H")

print(paste(c("mineral SOC", "organic SOC")[i],
            " at equilibrium =", sum(socE), sep = ""))
}

##time series
years=seq(from=1,to=100,by=1)
Litter=data.frame(year=years,Litter=rnorm(n=length(years),mean=5,sd=2))
ClimData=data.frame(years,Temp=7+sin(2*pi*years)+rnorm(n=length(years),mean=0,sd=1),
                    Moisture.Min.Welldrained=15+sin(2*pi*years)+rnorm(n=length(years),mean=0,sd=5),
                    Moisture.Org.Paludified=60+sin(2*pi*years)+rnorm(n=length(years),mean=0,sd=5))


print("Time series data (yearly):")
print(str(cbind(ClimData,Litter)))
print(c("AWEN=",AWEN))


#environmental effect
#this can be replaced by any soil TEMP, SWC function
y07.EMfun <- function(MT,SWC){ 
  qval=1.743
  wopt = 31.69
  EM= 0.999^((wopt - SWC)^2)*qval^((MT-10)/10)
  EM #environmental modifier
}

par(mfrow=c(1,2))
matplot(ClimData[,3:4],type="l", lty=1, col=1:2, xlab="Years", ylab="SWC min and organic")
plot(y07.EMfun(SoilTemperature,ClimData[,3]), type="l", ylim = c(0,1), 
     ylab="environmental modifier", xlab = "years")
lines(y07.EMfun(SoilTemperature,ClimData[,4]), col = 2)
legend("topleft",c("mineral","organic"), col = c(1,2), lty =1, bty="n")

##run the model and plot carbon pools and respiration
par(mfrow=c(2,2), mar=c(4,5,3,1))
for(i in 1:2){

  #i = 1
  SoilTemperature <- ClimData$Temp 
  Soilwatercontent<- ClimData[,c("Moisture.Min.Welldrained","Moisture.Org.Paludified")[i]]
  
  yassomio.ts <- Yasso07mio.EYM(t.step=1, #years
                                t=years, #years
                                C0=rep(0,5), #initial carbon
                                AWEN = c(0.52,0.18,0.08,0.2,0), #to separate litter to yasso AWEN pools, this depends on plant organ and species
                                In=Litter,#litter C input (same length as years)
                                xi = 0, # only xi = 1  will replace climate data no climate effect,
                                MT= SoilTemperature,
                                SWC=Soilwatercontent,
                                WS=0, #WS 2 = woody e.g. branches, coarse roots
                                ym.mean.params) 

Ct=yassomio.ts[["Ct"]]
Rt=yassomio.ts[["Rt"]] #respiration

#plot carbon pools
matplot(Ct, type="l", lty=1, col=1:5, ylim=c(0,max(rowSums(Ct))), 
        xlab="Years", ylab="Carbon stocks", main =c("Moisture.Min.Welldrained","Moisture.Org.Paludified")[i])
lines(rowSums(Ct), lty=1, col=1, lwd=2)
legend("topleft", c("A","W","E","N", "H"), lty=1, col=1:5, bty="n", n = 2)
#plot respiration
matplot(rowSums(Rt), type="l", ylab="Respiration", lty=1, col=1:2)

}

## Litter input in form of AWEN for months
months <- years
AWEN <- c(0.52,0.18,0.08,0.2,0)
Litter.months=data.frame(month=months,Litter.month=rnorm(n=length(months),mean=(5/12),sd=1/12)) #note yearly litter was 10, monthly 10/12
LA =  matrix(AWEN , nrow=length(months), ncol=5, byrow=TRUE) 
LI.awen = as.data.frame(LA*as.vector(Litter.months[,2])) #LitterInput
names(LI.awen) = c("A","W","E","N", "H")
LI.y.awen = data.frame(months, LI.awen)
as.matrix(LI.y.awen[,2:6])

#soil temp and swc have the same length 100, and are used from above
ClimData.months <- ClimData
names(ClimData.months)[1]<- "months"

print("Time series data (monthly):")
print(str(cbind(ClimData.months,Litter.months)))
print(c("AWEN=",AWEN))



yassomio.ts <- Yasso07mio.EYM(t.step=1/12, #months
                              t=months, #months
                              C0=rep(0,5), #initial carbon
                              AWEN = c(0,0,0,0,0), # or AWEN e.g. c(0.52,0.18,0.08,0.2,0), #to separate litter to yasso AWEN pools, this depends on plant organ and species
                              In=LI.y.awen ,#litter C input (same length as years)
                              xi = 0, # only xi = 1  will replace climate data no climate effect,
                              MT= SoilTemperature,
                              SWC=Soilwatercontent,
                              WS=0, #WS 2 = woody e.g. branches, coarse roots
                              ym.mean.params) 

Ct=yassomio.ts[["Ct"]]
Rt=yassomio.ts[["Rt"]] #respiration

#plot carbon pools
matplot(Ct, type="l", lty=1, col=1:5, ylim=c(0,max(rowSums(Ct))), 
        xlab="Months", ylab="Carbon stocks", main =c("Moisture.Min.Welldrained","Moisture.Org.Paludified")[i])
lines(rowSums(Ct), type="l", lty=1, col=1, lwd=2)
legend("topleft", c("A","W","E","N", "H"), lty=1, col=1:5, bty="n", n = 2)
#plot respiration
matplot(rowSums(Rt), type="l", ylab="Respiration", xlab="Months",lty=1, col=1:2)

}
cat("\n Yasso07mio model: the model function 'Yasso07mio.EYM' and parameters 'ym.mean.params' were loaded, for simple examples of data input and model simulation run Yasso07mio.EYM.examples() ")
