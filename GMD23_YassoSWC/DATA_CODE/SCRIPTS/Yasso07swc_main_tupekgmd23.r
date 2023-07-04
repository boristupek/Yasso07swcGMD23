
##################################################################################################################
# Code to replicate the analysis of Tupek et al. 2023: 
# Modeling boreal forest’s mineral soil and peat C stock dynamics with Yasso07 model coupled with updated moisture modifier. 
# Geoscientific model development
# 
# 
#
# I) Run Yasso07 soil C model in original form for mineral soils (Tuomi et al. 2011)
#                                              and for peatlands with wetlands reduction (Kleinen et al. 2021)
# II) Run Yasso07 model coupled with updated soil temp and moisture functions and 
#     calibrated with data from Vatiharju-Lakkasuo boreal frorest - mire ecotone in Finland (Tupek et al. 2008, Tupek et al. 2015)
###################################################################################################################

# Boris Tupek
# Natural Resources Institute Finland (LUKE)
# June 2023
# boris.tupek@luke.fi


# CODE STRUCTURE:
# 
# 1) plot Figure 2 (boxplots) and Figure 3 (time-series) data
#    of the daily CO2 measurements, soil temperature and moisture from 2004-2006 period,
#    and soil C stocks and litter data from 9 sites of forest-mire ecotone
#
# 2) run Yasso07 model in original version for mineral soils (Tuomi et al. 2011) and with wetland reduction (Kleinen et al. 2021)
#    based on litter input data (in total and separated to its AWEN quality fractions) (plot Figures S2 and S3) 
#             and meteorological data
#
# 3) MCMC Bayesian data assimilation
#    run Yasso07 model with modified environmental functions
#    plot Figure 4 (calibrated temp and swc modifiers)
#    plot Figure 5 (measured and modeled SOC and CO2, and their residual trends)
#
# 4) Evaluate the Yasso07 original model performance against the Yasso07 with updated t,swc modifier

rm(list=ls())

#note: modify paths to your own workspace! 
path.data <- "D:/LUKE/GMD23_YassoSWC/DATA_CODE/DATA/"
path.results <- "D:/LUKE/GMD23_YassoSWC/DATA_CODE/RESULTS/"
path.scripts <- "D:/LUKE/GMD23_YassoSWC/DATA_CODE/SCRIPTS/"


## preparation objects for plotting the figures#########
## color palettes
#install.packages('squash', dependencies = T)
library(squash)
#define colors
blackwhite.palette <- colorRampPalette(c("white", "black"), space = "rgb")
bw.colors <- blackwhite.palette(110)[10:110]
bw.palette <- colorRampPalette(c("black", "white"), space = "rgb")
bluegreen.palette <- colorRampPalette(c("blue", "green"), space = "rgb")
greenblue.palette <- colorRampPalette(c("green", "lightblue"), space = "rgb")
gb.colors <- greenblue.palette(110)[10:110]


# colors for 3 groups, forest, transitions, mires
ge.clrs <- colorRampPalette(c("white", "darkgreen"))
b.clrs <- colorRampPalette(c("white", "darkblue"))
g.clrs <- colorRampPalette(c("white", "black"))
ge5.clrs <- ge.clrs(12)[c(1,3,6,9,11)]
b5.clrs <- b.clrs(12)[c(1,3,6,9,11)]
g5.clrs <- g.clrs(12)[c(1,3,6,11,9)]
g5.eq.clrs <- g.clrs(12)[c(1,3,6,11,9,1)]

#A function to add arrows on the chart ###########
error.bar <- function(x, y, upper, lower, length=0.1,...){ #lower=upper, 
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

#############################################
## PLOTS OF MEASURED DATA #####################################################################
# Rh,T,SWC data ################################################################################
r <- read.delim(paste0(path.data,"r.notveg_t5sm10_eco.g.ft.csv"), header = T, skip = 0, sep=";")
names(r)
#observed soil respiration for 9 sites of the forest-mire ecotone
rno <- r[,c("ft","year","month","rn")]
rno<- rno[complete.cases(rno),]
names(rno)<-c("plot","year","month","rno")

# observed soil respiration WITH temp, swc data 
rno2 <- r[,c("ft","year","month","rn","t5mo1","sm10mo1")]
rno2<- rno2[complete.cases(rno2),]
names(rno2)<-c("plot","year","month","rno","t5mo1","sm10mo1")
#View(rno2)

## SOC interpolated for 1 m (FIG S1), and SOC of profiles for interpolation 
soilc.bg15 <- read.delim(paste0(path.data,"soilc.bg15.csv"), header = T, skip = 0, sep=";")
soilc.oct2015 <- read.delim(paste0(path.data,"soilc.oct2015.csv"), header = T, skip = 0, sep=";")
tSOC <- read.delim(paste0(path.data,"tSOC.csv"))$x
tSOCw <- read.delim(paste0(path.data,"tSOCw.csv"))$x
tSOC.ci1 <- read.delim(paste0(path.data,"tSOC.ci1.csv"))$x
tSOC.ci99<- read.delim(paste0(path.data,"tSOC.ci99.csv"))$x

SOC.e.obs <- data.frame(tSOC,tSOC.ci1,tSOC.ci99)
SOC.e.obs
#observed forest/mire socs
socobs <- tSOC
#monthly time series of observed soc (assumed insignificant SOC changes)
socobs.rep <- rep(socobs,each=36)

ls()

#plot Fig 2 ##
source(paste0(path.scripts,"y7swc_gmd23_fig2.r"))

# plot Fig 3 ##
source(paste0(path.scripts,"y7swc_gmd23_fig3.r"))

##################################################################################################################
## EMPIRICAL CO2 NON-LINEAR REGRESSION MODELLING #############################################################################################

#2004-2006 years models RH - soil heterotrophic respiration
print("###Parameters-statistics of TABLE 1 #####################################")
nls.pars9.list <- list()
#fit on forest/mire site level
for(i in 1:9){
  #i = 1
  r.i <- subset(r,ft == i)
  #mires
  nls.rn.i <- nls(rn ~ rref*0.998^((wopt - sm10mo1)^2)*q^((t5mo1-10)/10), 
                  data = r.i,
                  na.action= na.omit, start = c(rref= c(rep(0.48,4),rep(0.40,3),rep(0.31,2))[i],
                                                q = c(rep(2.65,4),rep(2.28,3),rep(1.37,2))[i],
                                                #NOTE!!! setting d param to 0.998 for easier convergence
                                                #d = 0.998,
                                                wopt = c(rep(15,4),rep(60,3),rep(60,2))[i]))
  
  nls.pars9.list[[i]] <- nls.rn.i
  print(c(i,c("CT","VT","MT","OMT","OMT+","KGK","KR","VSR1","VSR2")[i]))
  print(summary(nls.rn.i))
}
#fit on foret-mire ecotone level
nls.rn.eco <- nls(rn ~ rref*d^((wopt - sm10mo1)^2)*q^((t5mo1-10)/10), 
                data = r,
                na.action= na.omit, start = c(rref= 0.4,
                                              q = 2,
                                              d = 0.998,
                                              wopt = 30))
print("ECOTONE")
print(summary(nls.rn.eco))

##################################################################################################################
## SOC and CO2 PROCESS BASED MODELLING #############################################################################################
##################################################################################################################

#LITTER INPUT  (data to plot Fig S2 and Fig S3)###############################
#based on 
#stands <- read.csv(paste0(path.data,'Stands.ecotone.csv'), header=T, sep = ",",stringsAsFactors = FALSE, na.strings=c("NA",""))
#
#complete analysis of stand data to litter and AWEN fractions (RUN_LIT_AWEN = "YES" to run if needed)
RUN_LIT_AWEN = "NO" #LOAD ANALYSIS OUTPUTS  
if(RUN_LIT_AWEN == "YES"){ #RUN THE STAND TO LITTER-AWEN ANALYSIS
source(paste0(path.scripts,"biomass_litter_awenyasso_function.r"))
  #this requires 4 other Rscripts BIOMASS_functions.r,biomass.core.r,LITTER_functions.r,litter_updated.r 
  }else{ 
#read litter outputs of stand data analysis
# biomass and litter modeled by method of Finnish GHG described in Tupek et al 2019
load(file=paste0(path.data,"biomass.litter_trees.understory_years.RData"))
#includes litter "ul" understory,"tl"- trees, "tl.msd" agregated for species, "tul.sum" agregated for sites,
#"tul.bp" for plotting (NOTE: here tree finerots are excluded from the litter input as the tree roots were cut)

#load litter input from trenching for 2004 JULY Time series input!!!
load(paste0(path.data,"litter.trench.cut_coarse.fineroots.RData"))
#includes "tree.w.coarseroot.cut04.awen.sum", "treeunder.nw.fineroot.cut04.awen"

#load lists of simulated AWEN data woody and nowoody litter 
load(file=paste0(path.data,"litter.eco_equilibrium_months.RData"))
#includes lists "litter.eq.awen", "litter.natural.month.awen","litter.fc.month.awen"
  }

#woody and non-woody awen litter
treeunder.w.eq.awen <- litter.eq.awen[[1]]
treeunder.nw.eq.awen <- litter.eq.awen[[2]]

#total monthly litter input
litter.fc.month.awen_tot <-litter.fc.month.awen[[5]]

## plot Fig. S2 litter barplots
source(paste0(path.scripts,"y7swc_gmd23_fig_S2.r"))

## plot Fig. S3 litter timeseries
source(paste0(path.scripts,"y7swc_gmd23_fig_S3.r"))


#MICROMETEOROLOGICAL DATA ###########################################

## TEMP and WATER CONTINUOUS DAILY,MONTHLY,ANNUAL DATA ###
## Tair and Precipitation for Yasso07 original
load(file=paste0(path.data,"TA.PREC.TAMP.eco_equilibrium_ymd.RData"))
ls()[grep("ta.prec",ls())]
#[1] "ta.prec.doy.3456"        "ta.prec.doy.456"         "ta.prec.month.3456"     
#[4] "ta.prec.month.456"       "ta.prec.taamp.year.3456" "ta.prec.taamp.year.456" 
#[7] "ta.prec.year7403.mean"  

## Tsoil and SWC for Yasso07 modified
load(file=paste0(path.data,"T5.SWC10.eco9_dm.RData"))
#includes "tn5.h","t5.h","sm10.h","tn5.m","t5.m","sm10.m" #tn5.h and t5.h identical, h hour, m month
#View(tn5.h)
## agreagate site level data for the SOC+CO2 plot below 
swc.9.mean <- mapply(FUN = function(x) { c(mean(x, na.rm=TRUE))},sm10.h[,grep("swc10",names(sm10.h))])
swc.9.sd <- mapply(FUN = function(x) { c(sd(x, na.rm=TRUE))},sm10.h[,grep("swc10",names(sm10.h))])
swc.9.min <- mapply(FUN = function(x) { c(min(x, na.rm=TRUE))},sm10.h[,grep("swc10",names(sm10.h))])
swc.9.max <- mapply(FUN = function(x) { c(max(x, na.rm=TRUE))},sm10.h[,grep("swc10",names(sm10.h))])

## agreagate site level data for the SOC+CO2 plot below 
tn5.9.mean <- mapply(FUN = function(x) { c(mean(x, na.rm=TRUE))},tn5.h[,grep("tn5",names(tn5.h))])
tn5.9.sd <- mapply(FUN = function(x) { c(sd(x, na.rm=TRUE))},tn5.h[,grep("tn5",names(tn5.h))])
tn5.9.min <- mapply(FUN = function(x) { c(min(x, na.rm=TRUE))},tn5.h[,grep("tn5",names(tn5.h))])
tn5.9.max <- mapply(FUN = function(x) { c(max(x, na.rm=TRUE))},tn5.h[,grep("tn5",names(tn5.h))])
#standard error
se <- function(x) sqrt(var(x, na.rm = T)/length(x))
tn5.9.se <- mapply(FUN = function(x) { sqrt(var(x, na.rm = T)/length(x))},tn5.h[,grep("tn5",names(tn5.h))])
swc.9.se <- mapply(FUN = function(x) { sqrt(var(x, na.rm = T)/length(x))},sm10.h[,grep("swc10",names(sm10.h))])

clim.tn5.swc.e9 <- data.frame(tn5=tn5.9.mean,swc=swc.9.mean)

###############################################################################################
### YASSO07 ORIGINAL (TUOMI et al 2011) #######################################################

print("#### YASSO07 model function ##############################################################")
#source yasso07 model function (correct version with litter input to all A,W,E,N pools compatible with SoilR, not just into A which is erroneous)
source(paste0(path.scripts,"yasso07Model.soilr.fi_function.r"))
#this includes function for equilibrium

#yasso07 function for monthly time series 
source(paste(path.scripts,"yasso07Model.soilr.fi_month_function.r",sep=""))
# yasso.month.example() #execute for an example run

#list functions
lsf.str()

#define parameter ranges for yasso07
#from Tuomi et l. 2011 
y07.params.names<-c("kA","kW","kE","kN","aA_W","aA_E","aA_N","aW_A","aW_E","aW_N","aE_A","aE_W","aE_N","aN_A","aN_W","aN_E","b1","b2","G","pH","kH","d1","d2","r")
y07.mean.params <- c(0.73,	5.8,	0.29,	0.031,	0.48,	0.01,	0.83,	0.99,	0,	0.01,	0,	0,	0.03,	0,	0.01,	0.92,	0.096,	-1.4,	-1.21,	0.0045,	0.0017,	-1.71,	0.86,	-0.306)
y07.min.params <- c(0.62,	5,	0.24,	0.027,	0.41,	0,	0.6,	0.94,	0,	0,	0,	0,	0,	0,	0,	0.79,	0.078,	-2.4,	-1.06,	0.0037,	0.0014,	-1.9,	0.76, -0.321) 
y07.max.params <- c(0.84,	6.6,	0.35,	0.042,	0.54,	0.16,	0.98,	1,	0.08,	0.21,	0.004,	0.003,	0.25,	0.007,	0.031,	0.99,	0.122,	-0.8,	-1.36,	0.0056,	0.0019,	-1.5,	0.96,	-0.29)
y07.params <- data.frame(n=1:length(y07.params.names),
                         names=y07.params.names,
                         min=y07.min.params,
                         mean=y07.mean.params,
                         max=y07.max.params)
y07.params

### RUN ORIGINAL YASSO07 (TUOMI et al 2011) ###########################################
## modeled using  original Yasso07 with reduction constant for wetlands
## Equilibrium SOC from 1 POOL of woody and non-woody litter #

# climate equilibrium (30 years)
clim.e <- as.numeric(ta.prec.year7403.mean)

# litter input
#View(litter.eq.awen)
#names(litter.eq.awen)
#[1] "treeunder.nw.eq.awen" "treeunder.w.eq.awen" 
treeunder.w.eq.awen <- litter.eq.awen[[1]]
treeunder.nw.eq.awen <- litter.eq.awen[[2]]

# EQUILIBRIUM Yasso simulations ############ 
#
# note: 
# Yasso07 runs as Tuomi (2011) for mineral, and Kleinen (2021) wetlands reduction for VSR1 and VSR2 !!!
# wetlands="y" YES for 8,9  sets wetlands reduction 

soc.awenh.e <- data.frame(matrix(NA,9,5))
names(soc.awenh.e) <- c("A", "W", "E","N", "H")
socmo1.e <- rep(NA,1)
rh.molitt1.e <- rep(NA,1)
for(i in 1:9){
  #i = 1
  inawen.wie  <- c(as.numeric(treeunder.w.eq.awen[i,2:5]),0)/10000 ##convert from ha to kg m2!!!!
  inawen.nwie <- c(as.numeric(treeunder.nw.eq.awen[i,2:5]),0)/10000 ##convert from ha to kg m2!!!!
  inawen.tie <- inawen.wie+inawen.nwie #total
  c.size <- 0 
  #yasso07 equilibrium ~ c.size 
  #note: model matrix depends on litter input via the size parameter which could be the proportion of 2*w/(w+nw)?
  #however, SOC equilibrium from separated W and NW litter above, and pooled together were comparable ONLY if c.size = 0, see TEMP.soc.co2.mod.bayes.r
  if(i<8){
    AYS.ix <- yasso.matrix.fun(WS =c.size, clim =  clim.e, wetlands ="n", A.print = "n") 
  }else{
    AYS.ix <- yasso.matrix.fun(WS =c.size, clim =  clim.e, wetlands ="y", A.print = "n") 
  }
  soileqc.ix = -1*solve(AYS.ix) %*% inawen.tie # inawen.tie = TOTAL litter #inverse of matrix solve(B)
  soileqc.ix 
  
  soc.awenh.e[i,] <- soileqc.ix 
  socmo1.e[i] <- sum(soileqc.ix)
  rh.molitt1.e[i] <- sum(inawen.tie)
}
#rename outputs of orig model with wetlands reduction as "TW"
soc.awenh.TW.e <-soc.awenh.e
socmo1.TW.e <- socmo1.e
rh.molitt1.TW.e <- rh.molitt1.e

# MONTHLY Yasso simulations ############
# note: for 8,9 set for wetlands reduction wetlands="y" 
clim.m.e <- ta.prec.month.3456 #monthly weather in Juupajoki during 2004:2006
dim(clim.m.e)
for(i in 1:9){
  #i = 8
  litt.fcm.i <- subset(litter.fc.month.awen_tot, plot== i)
  #ix.jl04 <-which(litt.fcm.i$year==2004 & litt.fcm.i$month==7)
  #sum(litt.fcm.i[ix.jl04, c("A.tot", "W.tot", "E.tot", "N.tot")])/10000
  #2.49788 #note e.g. in site 1 adding 2.5 kg by trenching, harvesting
  litt.fcm.awen.i <- litt.fcm.i[,c("A.tot", "W.tot", "E.tot", "N.tot")] #litter on montly level!!!
  
  #compare model spinup RH to  sum of litter
  C0.i = as.numeric(soc.awenh.e[i,])
  n.i = dim(litt.fcm.i)[1]
  months.i = 1:n.i 
  MT = clim.m.e[,"tair.mean"] # should be of n.i length
  PR_mm = clim.m.e[,"precip.sum"]
  #Note: no temp amplitude for monthly temperatures, as they are used directly in modif temp function
  c.size <- 0
  
  inawen.ie <- data.frame(matrix(NA,n.i,6))
  names(inawen.ie) <- c("months","A", "W", "E","N", "H")
  inawen.ie$months <- months.i
  inawen.ie[,c("A", "W", "E","N")] <- litt.fcm.awen.i[,c( "A.tot","W.tot","E.tot","N.tot")]/10000 #convert from ha to kg m2!!!! 
  inawen.ie$H <- 0
  LI <-as.matrix(inawen.ie)
  
  if(i < 8){ #FT 1:7 mineral soils, wetlands = "n"
    #note yasso decomposition rates k/12
    mod.soilc.months.ix <-Yasso07Modelfi.month(months.i, #1:n of years 
                                               C0=C0.i, #initial carbon #vector of five
                                               AWEN = 0, #0 if litter in is in awen form
                                               In=LI, #litter C input in AWEN form (same length as years)
                                               xi = 0, #0 to pecify climate, 1 ignore climate
                                               MT=MT,#MeanTemperature (same length as years)
                                               PR_mm=PR_mm,#Precipitation_mm (same length as years)
                                               wetlands = "n", # "n"  for NO 
                                               WS=0) #scalar 2 woody or 0 nonwoody, 0 for total litter
  }else{  #FT 8:9 wetlands, wetlands = "n" or "y"
    mod.soilc.months.ix <-Yasso07Modelfi.month(months.i, 
                                               C0=C0.i, 
                                               AWEN = 0,
                                               In=LI, 
                                               xi = 0, 
                                               MT=MT,
                                               PR_mm=PR_mm,
                                               wetlands = "y", # "y" for YES
                                               WS=0) 
  }
  soc.months.ix = data.frame(getC(mod.soilc.months.ix))
  rh.months.ix = data.frame(getReleaseFlux(mod.soilc.months.ix))#respiration
  names(soc.months.ix) <-  c("A", "W", "E","N", "H")
  names(rh.months.ix) <-  c("A", "W", "E","N", "H") 
  
  soc.months.i <- cbind(litt.fcm.i[,c("plot", "year", "month")],soc.months.ix)
  rh.months.i <- cbind(litt.fcm.i[,c("plot", "year", "month")],rh.months.ix)
  
  if(i == 1){
    soc.months.e9 <- soc.months.i
    rh.months.e9 <- rh.months.i
  }else{
    soc.months.e9 <- rbind(soc.months.e9,soc.months.i)
    rh.months.e9 <- rbind(rh.months.e9,rh.months.i)
  }
}

#2) Run after following 1st yasso month loop
#yasso months runs FT1:7 mineral, 8:9 wetlands (wetlands="y")
soc.months.TW.e9 <- soc.months.e9 
rh.months.TW.e9 <- rh.months.e9

soc.months.TW.e9 #rowSums(soc.months.TW.e9[,c("A","W","E","N", "H")]  #soc monthly time series
rh.months.TW.e9 #soc rh monthly time series

#yasso.orig.wet (MIN soils and WETLANDS reductions)
soc.twe <- socmo1.TW.e #modeled 
socres.twe<- socobs - soc.twe #residuals
socnr.twe <- socres.twe/socobs #n

#YASSO modeled Rh with TW  version
#(Tuomi 2011 for mineral, Kleinen 2021 wetlands reduction for VSR1 and VSR2 )!!!
#modeled respiration
rh.months.TW.e9$rh.m.tot<- rowSums(rh.months.TW.e9[, c("A", "W", "E", "N")])*44/12 #convert from Kg C m-2 month-1 Kg CO2 m2 month-1
rh.months.TW.e9.0 <- subset(rh.months.TW.e9, year>2003) #exclude 2003 for clarity in plots
rh.mm1 <-rh.months.TW.e9.0[,c("year","month", "plot", "rh.m.tot")]
rh.m.tot <- rh.months.TW.e9.0$rh.m.tot

#View(rh.mm1) #modeled
#View(rno2) #measured

#merge modeled rh (m -modeled) with observed: rn-respiration not vegetated o-observed 
rnom.twe <- merge(rno2, rh.mm1, by =c("plot","year","month"), all.x=T ) 
# View(rnom.twe) 
names(rnom.twe)
rnom.twe$rhmo_h <- rnom.twe$rh.m.tot*1000/(30*24) #convert modeled respiration from kg/monthly to g/hourly values
rnom.twe$rh_h.res <- rnom.twe$rhmo_h - rnom.twe$rno 
rnom.twe$norm.rh_h.res <- rnom.twe$rhmo_h/rnom.twe$rno-1 

co2h.obs <- rnom.twe$rno
co2h.twe <- rnom.twe$rhmo_h
co2h.twe.res <- rnom.twe$rh_h.res
co2h.twe.resn <- rnom.twe$norm.rh_h.res 
t5h.twe <-rnom.twe$tn5
swch.twe <-rnom.twe$swc

##################################################################################################################
### YASSO07 MODIFIED ENVIRONMENTAL EFFECT for mineral and organic soils (MIO) #####################################
print("#### TABLE 2 # MCMC calibration statistics ##############################################################")

#load yasso07mio.EMY function with new soil temp and moisture modifier 
#     for SOC and CO2 of equilibrium, and yearly/monthly time series 
source(paste0(path.scripts,"yasso07mio_Eym_function.r"))

#yasso07 function with new modifier running with MONTHLY LITTER INPUT AND MONTHLY DECAY RATES
## two equilibrium functions IDENTICAL = NOT SENSITIVE FOR THE TIME STEP
source(paste0(path.scripts,"yasso07mio_Eym_function_month.r"))

#list functions
lsf.str()

# parameters if Tuomi et al. (2011), except for environmental climate parameters
## climate effect parameters pE(d=0.999, q=1.743,wopt = 31.69)
#read mean params with lower and upper bounds
ym.params.lbu <- read.delim(paste0(path.scripts,"Yasso07params_ts_swc.csv"), header = T, skip = 0, sep=",")
ym.params.lbu

refPars.upper <- ym.params.lbu[,"upper"]
refPars.upper[23:24] <- c(5,70) #replace upper limits for q and wopt
refPars.best <- ym.params.lbu[,"best"]
refPars.lower <- ym.params.lbu[,"lower"]
refPars.blu <- data.frame(best= refPars.best, lower= refPars.lower, upper =refPars.upper)
er.mcmc <-data.frame(best=c(80,1,0.15,0.5),
                     lower=rep(0.0001, 4),
                     upper=c(150,3,0.7,3)) 
row.names(er.mcmc)<-c("soc.a","soc.b", "co2.a", "co2.b") # #parameters for the error terms
refPars.blu0 <- rbind(refPars.blu,er.mcmc) 
refPars.blu.soc <-refPars.blu0[1:26,]


## BAYESIAN MCMC of Yasoo07 with t,swc modifier  ###############################
library(SoilR) 
library(sensitivity)
library(BayesianTools)
library(DHARMa)
library(rgl)

## Load calibration functions for SOC and SOC-CO2 ###
source(paste0(path.scripts,"yassomioE.calibr.fun.D12_SOC_SOCCO2h.r"))

#selection of calibrated parameters (d,q10,wopt + 2errs for soc or 4errs for socco2)
parind.soc <- 22:26 #SOC
parind0 <- 22:28 #SOCCO2

# SOC likelihood ###############################
likelihood2.soc <- function(parameter, sum = T){
  x = refPars.blu.soc$best
  x[parind.soc] <-  parameter
  
  ymio.out.SOC.x <- yassomioE.calibr.fun.D12.SOC(x)
  soc.e.mod <- ymio.out.SOC.x[["socmio1.e"]]
  
  #residuals
  socres <- socobs-soc.e.mod
  
  llvalues.soc <- sum(dexp(abs(socres), 
                           rate = x[25]+x[26]*soc.e.mod,log=T)) 
  return(llvalues.soc)
}
#test SOC likelihood
refPars.best <- refPars.blu.soc$best
true = refPars.best[parind.soc]
true
#test 
likelihood2.soc(true)
likelihood2.soc(true*1.0001)

# RUN BAYESIAN SOC calibration #####
BAYES_SOC = "NO" # skip by default and load saved chains 
if (BAYES_SOC == "YES"){ #if YES on an average laptop pc runs this for approx 50 min with MCMC iterations set to 30000    
# Define prior
prior.soc <- createUniformPrior(lower =  refPars.blu.soc$lower[parind.soc], 
                                upper = refPars.blu.soc$upper[parind.soc])
# Bayes test setup 
BSmod1.soc <- createBayesianSetup(likelihood2.soc, prior.soc, 
                                  best = refPars.blu.soc$best[parind.soc], 
                                  names = rownames(refPars.blu.soc)[parind.soc],
                                  parallel = 9) #9 cores, check pc system if 9 is OK
summary(BSmod1.soc)
# Define running the sampler
# for test run (model will not converge) reduce number of iterations to e.g. 1000 = approx 3 min 
settings1 <- list(iterations = 30000,
                  nrChains = 3, message = TRUE)
#run the MCMC
chain1.soc <- runMCMC(BSmod1.soc, sampler="DEzs", settings = settings1)
# save the 'dismo::ModelEvalution' object
saveRDS(chain1.soc, file=paste0(path.results,"bayes.fit.SOC.30000.3chains.rds"))	
#diagnostics
gelmanDiagnostics(chain1.soc) $mpsrf #OK if around 1
#[1] 1.059662 
summary(chain1.soc)
# MCMC sampler:  DEzs 
# Nr. Chains:  9 
# Iterations per chain:  10000 
# Rejection rate:  0.926 
# Effective sample size:  1193 
# Runtime:  2603.73  sec. 

# Parameters
#psf   MAP    2.5% median  97.5%
#22    1.027 0.999 0.999  0.999  1.000
#23    1.033 4.364 0.795  3.959  4.951
#24    1.028 6.762 5.182  8.416 60.703
#soc.a 1.029 0.145 0.003  0.098  1.584
#soc.b 1.026 0.000 0.000  0.000  0.018

## DIC:  341.408 
## Convergence 
#  Gelman Rubin multivariate psrf:  1.06 
plot(chain1.soc)

}
# end of bayes soc calibration ###########

# SOC and hourly CO2 (Rh) likelihood #################################
likelihood3.SOC.CO2hour <- function(parameter, sum = T){
  x = refPars.blu0$best
  x[parind0] <- parameter 
  #monthly
  ymio.out.SOCCO2.x <- yassomioM.calibr.fun.D12.SOC.CO2month(x)
  
  soc.month.mod <- ymio.out.SOCCO2.x[["soc.mio.months.e9"]]$socM.i
  co2.month.mod <- ymio.out.SOCCO2.x[["rh.mio.months.e9"]]$rhM.i
  
  #merge modeled rh (m -modeled) with observed: r-respiration n-not_vegetated o-observed 
  rnom.x <- merge(rno, ymio.out.SOCCO2.x[["rh.mio.months.e9"]], by =c("plot","year","month"), all.x=T )
  rnom.x$rhM.i_h <- rnom.x$rhM.i*1000/(30*24) #convert modeled respiration from kg/monthly to g/hourly values
  
  ##residuals (mod-obs)
  socres_m <- (soc.month.mod - socobs.rep) 
  rnom.x$rno[which(rnom.x$rno==0)]<-NA
  co2res_h <- rnom.x$rhM.i_h-rnom.x$rno
  
  #SOC 
  llvalues.soc.month <- sum(dexp(abs(socres_m), rate=1/(x[25]+x[26]*soc.month.mod),log=T),na.rm=T)
  #Rh
  llvalues.co2.hour <-sum(dexp(abs(co2res_h), rate = 1/(x[27]+x[28]*rnom.x$rhM.i_h),log=T), na.rm=T)
  #CO2 10 % weight
  llvalues = llvalues.soc.month + 0.1*llvalues.co2.hour 
  
  return(llvalues)
}
# SOC-CO2 calibration # same as for SOC but replace the SOC likelihood by SOCCO2 (note computation takes about 5 hours)

#MCMC sampling chains, the runs here are loaded to save the computation time
# SOC fit as in soc.co2.mod.bayes.brt_17.03.23.r
# #load the sampling likelihood2.soc
chain1.x.soc <- readRDS(paste0(path.results,"bayes.fit.SOC_lin.errpars.30000.3chains.eco_15.03.23.rds")) 

# SOCCO2 fit as in soc.co2.mod.bayes.brt_17.03.23.r
# #load the sampling likelihood3.SOC.CO2hour weight 10% CO2
chain1.x0 <- readRDS(paste0(path.results,"bayes.fit.SOCCO2hour.varmod.errpars.like3.30000.3chains.eco_17.03.23.rds"))

## chain statistics #
summary(chain1.x.soc) #SOC
# plot(chain1.x.soc)
summary(chain1.x0) #SOCCO2
# plot(chain1.x0)

# extract MAP (mean aposteriori parameters) from the MCMC chains 
map.soc <- MAP(chain1.x.soc)$parametersMAP[c("d","q","wopt")] 
map.socco2 <- MAP(chain1.x0)$parametersMAP[c("d","q","wopt")]

## plot Fig. 4 environmental modifiers with MAP ################################
source(paste0(path.scripts,"y7swc_gmd23_fig4.r"))


### RUN CALIBRATED YASSO07 with T,SWC function ###############################################
#replace q d wotp params by MAP from bayesian MCMC calibration on ecotone level ##############
x.pars = ym.mean.params 
for(i in 1:2){
  # i = 1
  x.pars[c(22:24)] <- as.numeric(get(c("map.soc", "map.socco2")[i])) 
  
  #SOC and CO2 simulations
  ymio.out.SOCCO2.x <- yassomioM.calibr.fun.D12.SOC.CO2month(x.pars)
  soc.month.mod0 <- ymio.out.SOCCO2.x[["soc.mio.months.e9"]]
  #equilibrium SOC
  soc.mio <- subset(soc.month.mod0, year==2004 & month==1)$socM.i
  #monthly SOC and CO2
  soc.month.mod <- ymio.out.SOCCO2.x[["soc.mio.months.e9"]]$socM.i
  co2.month.mod <- ymio.out.SOCCO2.x[["rh.mio.months.e9"]]$rhM.i
  
  #merge modeled rh (m -modeled) with observed: rn-respiration not vegetated o-observed 
  rnom.x <- merge(rno2, ymio.out.SOCCO2.x[["rh.mio.months.e9"]], by =c("plot","year","month"), all.x=T )
  rnom.x$rhM.i_h <- rnom.x$rhM.i*1000/(30*24) #convert modeled respiration from kg/monthly to g/hourly values
  
  ## add CO2 residuals (absolute and relative) 
  rnom.x$rno[which(rnom.x$rno==0)]<-NA
  rnom.x$co2res_h <- rnom.x$rhM.i_h-rnom.x$rno #absolute
  rnom.x$norm.co2res_h <- rnom.x$rhM.i_h/rnom.x$rno-1 #relative
  
  #SOC and CO2 model outputs based on parameter fitting of 1 and 2
  #1 p(theta,SOC) , 2 p(theta,SOCCO2)
  assign(c("soc.de1","soc.de2")[i],soc.mio) 
  assign(c("memo.co2h.de1","memo.co2h.de2")[i],rnom.x) 
  assign(c("co2h.de1","co2h.de2")[i],rnom.x$rhM.i_h) 
}

names(memo.co2h.de1)

#residuals
#1) #soc ecotone ~ SOC
socres.de1 <- socobs - soc.de1 #residuals
socnr.de1 <- socres.de1/socobs #normalized residuals
#co2 ecotone  ~ SOC
co2res.de1 <- memo.co2h.de1$co2res_h #co2obs - co2.de1 
co2nr.de1 <- memo.co2h.de1$norm.co2res_h#co2res.de1/co2obs 

#2) #soc ecotone ~ SOCCO2
socres.de2 <- socobs - soc.de2 #residuals
socnr.de2 <- socres.de2/socobs #normalized residuals
#co2 ecotone  ~ SOCCO2
co2res.de2 <- memo.co2h.de2$co2res_h #co2obs - co2.de1 
co2nr.de2 <- memo.co2h.de2$norm.co2res_h#co2res.de1/co2obs 

tn5.hour <- memo.co2h.de1$t5mo1 
swc.hour <- memo.co2h.de1$sm10mo1  

## plot Fig. 5 Scatterplots SOC and CO2 - measured vs modeled and residuals #####
source(paste0(path.scripts,"y7swc_gmd23_fig5.r"))

###########################################################################################################
## YASSO07 MODELS PERFORMANCE statistics  #################################################################

#parameters Yasso07 
npars.ytw = length(y07.params.names) #24 in Tuomi et al. 2011
npars.yd = 24 #in modified version here
npars.i <- 24 #set to 24 as the original and modified Yasso07 versions have same number of parameters

#To compare models using AIC, you need to calculate the AIC of each model. 
#If a model is more than 2 AIC units lower than another, then it is considered significantly better than that model.
mod.err <- data.frame(data=rep(c("SOC","CO2"),each=3),
                      model=rep(c("Y07.TW","Y07.Dsoc","Y07.Dsocco2"),2),
                      adj.r2=NA,
                      rmse_kg = NA,
                      aic = NA,
                      mbe = NA,
                      mae = NA)
for(j in 1:2){  
  #j=1
  obs <-get(c("socobs","co2h.obs")[j]) #use hourly CO2 
  for(i in 1:3){
    #i=1
    if(j==1){
      mod.y07.i <- get(c("soc.twe","soc.de1","soc.de2")[i])
    }else{
      mod.y07.i <- get(c("co2h.twe","co2h.de1","co2h.de2")[i]) 
    }
    fit.i <- summary(lm(obs ~ mod.y07.i))
    adj.r2.i <- fit.i[["adj.r.squared"]]
    mod.bias.i <- mod.y07.i - obs
    rmse_kg.i <- sqrt(sum((obs - mod.y07.i)^2, na.rm=T)/length(obs))
    aic.i <- length(obs)*log(sqrt(sum((obs - mod.y07.i)^2, na.rm=T)/length(obs))) + 2*npars.i
    mbe.i <- mean(mod.bias.i, na.rm=T)
    mae.i <- mean(abs(mod.bias.i), na.rm=T)
    if(j==1){
      mod.err[i,3:7]<-round(c(adj.r2.i,rmse_kg.i,aic.i,mbe.i,mae.i),2)
    }else{
      mod.err[(i+3),3:7]<-round(c(adj.r2.i,rmse_kg.i,aic.i,mbe.i,mae.i),2)
    }
  }
}
length(co2h.obs)
print("## TABLE 2 #############################################################")
print("Statistics of the Yasso07 original model performance against the Yasso07 with updated t,swc modifier")
print(mod.err)
