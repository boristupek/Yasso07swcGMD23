
##################################################################################################################
# Reconstruct figures of Tupek et al. 2023: 
# Modeling boreal forest’s mineral soil and peat C stock dynamics with Yasso07 model coupled with updated moisture modifier
#
###################################################################################################################

# Boris Tupek
# Natural Resources Institute Finland (LUKE)
# June 2023
# boris.tupek@luke.fi

#note: standalone requires r data and SOC data
# Rh,T,SWC data #################
#r <- read.delim(paste0(path.data,"r.notveg_t5sm10_eco.g.ft_19.06.23.csv", header = T, skip = 0, sep=";"))

# SOC data
#soilc.bg15 <- read.delim( paste0(path.data,"soilc.bg15_19.06.23.csv", header = T, skip = 0, sep=";"))
#soilc.oct2015 <- read.delim( paste0(path.data,"soilc.oct2015_19.06.23.csv", header = T, skip = 0, sep=";"))
#tSOC <- read.delim( paste0(path.data,"tSOC_19.06.23.csv", header = T, skip = 0, sep=";"))
#tSOCw <- read.delim(paste0(path.data,"tSOCw_19.06.23.csv", header = T, skip = 0, sep=";"))
#tSOC.ci1 <- read.delim( paste0(path.data,"tSOC.ci1_19.06.23.csv", header = T, skip = 0, sep=";"))
#tSOC.ci99<- read.delim(paste0(path.data,"tSOC.ci99.06.23.csv", header = T, skip = 0, sep=";"))


#################################################################################
## FIG 2 # BOXPLOTS OF SOC,SWC, RH ##############################################
## MEASURED SOC SOIL ORGANIC CARBON STOCKS ######################################
# combined data from Tupek et al, 2015 and October 2015 field campaign with Alla Yurova

#forest/mire SWC10 statistics 
swc.9.med <- boxplot(sm10mo1 ~ ft, data= r, plot=F)$stats[3,]
swc.9.min <- boxplot(sm10mo1 ~ ft, data= r, plot=F)$stats[1,]
swc.9.max <- boxplot(sm10mo1 ~ ft, data= r, plot=F)$stats[5,]


## SOC interpolated for 1 m (FIG S1), and SOC of profiles for interpolation 
#load(paste0(data.path,"SOC.ci_kgm2_ft9.RData"))

ls()
# SOC sums ~ SWC, CO2 - MS. FIG 2 Figure 
#figure.1m.SOC.ci_SWC.sd_CO2.boxplots()

par(mfcol=c(2,1), mar=c(0,0,0,0), oma=c(5,5,2,5))
#plot SOCS
plot(c(0.5,10.5),c(0,170), #ylim = rev(c(0,1000)),
     xaxt="n",
     cex.axis=1.3,  col = "white", xlab = "SOC (g/cm3)", ylab = "soil depth (cm)")
axis(1, at=seq(0.7,10.5,1.2),labels=NA, col.axis = 2)
barplot(tSOC, yaxt="n", add=T)
points(seq(0.7,10.5,1.2),tSOC,col=1,pch=16,cex=1.2)
arrows(seq(0.7,10.5,1.2), tSOC.ci1, seq(0.7,10.5,1.2), tSOC.ci99,
       length=0.05, angle=90, code=3, lty = 1, col = 1)
mtext(expression("SOC (kg C m"^{-2} ~ ")"), 2, 2.6, cex = 1.3)#Sigma ~"SOC
abline(v=c(4.9,8.5), col = 2, lty = 3)

#plot soil moisture on 2nd y axis
par(new = T)
par(col.lab= "blue")
plot(seq(0.7,10.5,1.2),rep(NA,9), type="p", col="white", pch=15, 
     yaxt ='n',xaxt ='n', 
     #xlim= c(0,1.3),
     ylim=c(0, 100),
     axes=F, xlab=NA, ylab=NA, cex=1)
points(seq(0.7,10.5,1.2),swc.9.med,col="blue",pch=15,cex=1.7)
arrows(seq(0.7,10.5,1.2), swc.9.min, seq(0.7,10.5,1.2), swc.9.max,
       length=0.05, angle=90, code=3, lwd=1,lty = 1, col = "blue")
axis(4, at=seq(0,100,20),labels=T,cex.axis=1.2, col = 4, col.axis= "blue")
mtext(expression(bar(SWC) ~"(%)"), 4, 3.1, cex = 1.3, col = 4)

#plot CO2
plot(c(0.5,10.5),c(0,2.1), #ylim = rev(c(0,1000)),
     xaxt="n", yaxt="n",cex.axis=1.3,  col = "white", xlab = "SOC (g/cm3)", ylab = "soil depth (cm)")
mtext(expression("R"[h]~~group("(",g~CO[2]~~m^{-2}~~hour^{-1},")")), 2, 2.5, cex = 1.3)
axis(1, at=seq(0.7,10.5,1.2),labels=c("CT","VT","MT","OMT","OMT+","KGK","KR","VSR1","VSR2"), 
     cex.axis=1.3)
axis(2, at=seq(0,2,0.5),labels=T, cex.axis = 1.3)
boxplot(rn ~ ft, data= r, xaxt="n", yaxt="n", at = seq(0.7,10.5,1.2), add=T)
mtext(expression("Forest/mire types"), 1, 3.3, cex = 1.3)
abline(v=c(4.9,8.5), col = 2, lty = 3)

mtext(expression("upland forests", "transitions", "mires"), 3, -3.3, col = 2, 
      outer = FALSE, at = c(2,7,9), cex = 1.3)

##################################################################################################################
