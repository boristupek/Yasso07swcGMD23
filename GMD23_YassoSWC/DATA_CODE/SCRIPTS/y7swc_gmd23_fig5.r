
##################################################################################################################
# Reconstruct figures of Tupek et al. 2023: 
# Modeling boreal forest’s mineral soil and peat C stock dynamics with Yasso07 model coupled with updated moisture modifier
#
###################################################################################################################

# Boris Tupek
# Natural Resources Institute Finland (LUKE)
# June 2023
# boris.tupek@luke.fi

#note: standalone requires r data loaded in "y7swc_gmd23.r" 
#
## needs runs of YASSO SOC and CO2 for original and modified-calibrated version 


################################################################################
## Plot FIG 5 SOC measured modeled and residuals ###############################
#
par(xpd=FALSE) #Keep R from graphing a line outside of a chart's area
ci.lines.fun <- function(x,y,datavect,col.line){
  model <- lm(y ~ x ) 
  confint(model, level=0.95)
  predicted.intervals <- predict(model,data.frame(x=datavect),interval='confidence',
                                 level=0.99)
  x1=datavect
  lines(x1,predicted.intervals[,1],col=col.line,lwd=2)
  lines(x1,predicted.intervals[,2],col=col.line,lwd=2, lty= 2)
  lines(x1,predicted.intervals[,3],col=col.line,lwd=2, lty= 2)
}

par(mfrow=c(2,3), mar= c(5,5,2,2), oma=c(1,1,1,1))
#SCATTERPLOT SOC measured and modeled
plot(socmo1.TW.e,tSOC,col="gray",pch=15, cex=1.5,xlim = c(0,170),ylim = c(0,170),  cex.axis=1.3, cex.lab=1.3,
     xlab=expression("SOC Yasso07 (kg C m"^{-2}~")"), ylab =expression("SOC Observed (kg C m"^{-2}~")"))
arrows(socmo1.TW.e, tSOC.ci1, socmo1.TW.e, tSOC.ci99,
       length=0.05, angle=90, code=3, lty = 1, col = "gray")
#points(socmo1.T.e,tSOC,col="gray",pch=0,cex=1.1)
#p(theta|SOC)
points(soc.de1,tSOC,col="orange",pch=16,cex=1.1)
arrows(soc.de1, tSOC.ci1, soc.de1, tSOC.ci99,
       length=0.05, angle=90, code=3, lty = 1, col = "orange")
#p(theta|SOC-CO2)
points(soc.de2,tSOC,col=1,pch=1,cex=1.1)
arrows(soc.de2, tSOC.ci1, soc.de2, tSOC.ci99,
       length=0.05, angle=90, code=3, lty = 1, col = 1)
abline(0,1, col =2, lty =2)
legend("bottom", legend=c(expression(Y07.~xi~TW), 
                          expression(Y07.~xi~D~"p("~theta~"|SOC)"),
                          expression(Y07.~xi~D~"p("~theta~"|SOC-CO2)")),
       pch = c(15,16,1), col = c("gray","orange",1), bg = "white", bty = "n",cex=1.3,  pt.cex = c(1.3,1.3,1))
legend("topright", "a)", bty = "n", cex = 1.5)

#plot SOC residuals vs TEMP
plot(tn5.9.mean, socnr.twe,  ylim = c(-1,1),   
     xaxt="n",# yaxt="n",
     cex.axis=1.3, cex.lab=1.3, pch=15, col = "gray",
     xlab = expression("T"[5]~"(°C)"),
     ylab = "norm. SOC residuals")
axis(1, at=seq(0,10,0.2),labels=T, cex.axis=1.3, cex.lab=1.3, col = 1)
#add mean trend and confidence lines
y <-socnr.twe
x <- tn5.9.mean
datavect <- seq(0,10,0.2)
ci.lines.fun(x,y,datavect,col.line="gray")
#
y <-socnr.de1
points(x, y, pch=16,  col ="orange")
ci.lines.fun(x,y,datavect,col.line="orange")
#
y <-socnr.de2
points(x, y, pch=1,  col =1)
ci.lines.fun(x,y,datavect,col.line=1)
legend("topright", "b)", bty = "n", cex = 1.5)

#plot SOC residuals vs SWC
plot(swc.9.mean, socnr.twe,  ylim = c(-1,1),   
     xaxt="n",
     cex.axis=1.3, cex.lab=1.3, pch=15, col = "gray",
     xlab = expression("SWC"[10]~"(%)"),
     ylab = "norm. SOC residuals")
axis(1, at=seq(0,100,20),labels=T, cex.axis=1.3, cex.lab=1.3, col = 1)#, col.axis= "blue")
#add mean trend and confidence lines
x <- swc.9.mean
datavect <- 1:100
y <-socnr.twe
ci.lines.fun(x,y,datavect,col.line="gray")
#
y <-socnr.de1
points(x, y, pch=16,  col ="orange")
ci.lines.fun(x,y,datavect,col.line="orange")
#
y <-socnr.de2
points(x, y, pch=1,  col =1)
ci.lines.fun(x,y,datavect,col.line=1)
legend("topright", "c)", bty = "n", cex = 1.5)
legend("bottomleft", legend=c(expression(Y07.~xi~TW), 
                              expression(Y07.~xi~D~"p("~theta~"|SOC)"),
                              expression(Y07.~xi~D~"p("~theta~"|SOC-CO2)")),
       pch = c(15,16,1), lty= 1, col = c("gray","orange",1), bg = "white", bty = "n",cex=1.3,  pt.cex = c(1.3,1.3,1))

library(scales)

#SCATTERPLOT CO2 measured vs modeledl
plot(co2h.twe, co2h.obs, col = alpha("gray", 0.3), #"gray",
     pch=15,cex=1.3, cex.lab=1.3,
     xlim = c(0,2),ylim = c(0,2),
     xlab=expression(Rh~Yasso07~~group("(",g~CO[2]~m^{-2}~hour^{-1},")")),
     ylab =expression(Rh~Observed~~group("(",g~CO[2]~m^{-2}~hour^{-1},")")) )
points(co2h.de1, co2h.obs, col = alpha("orange", 0.3),#col ="orange",
       pch=16)
points(co2h.de2, co2h.obs, col = alpha("black", 0.3), #col =1,
       pch=1)
abline(0,1, col =2, lty =2, lwd=1)
legend("topright", "d)", bty = "n", cex = 1.5)
legend("topleft", legend=c(expression(Y07.~xi~TW), 
                           expression(Y07.~xi~D~"p("~theta~"|SOC)"),
                           expression(Y07.~xi~D~"p("~theta~"|SOC-CO2)")),
       pch = c(15,16,1), col = c("gray","orange",1), bg = "white", bty = "n",cex=1.3,  pt.cex = c(1.3,1.3,1))

#plot normalized CO2 residuals vs TEMP
plot(tn5.hour, co2h.twe.resn,  ylim = c(-1,7),   
     xaxt="n",# yaxt="n",
     cex.axis=1.3, cex.lab=1.3, pch=15,col = alpha("gray", 0.3),
     xlab = expression("T"[5]~"(°C)"),
     ylab = expression(norm.~CO[2]~residuals))
axis(1, at=seq(0,20,5),labels=T, cex.axis=1.2, col = 1)
#add mean trend and confidence lines
x <- tn5.hour
x[is.na(x) | x=="Inf"] = NA
y0 <-co2h.twe.resn
y0[is.na(y0) | y0=="Inf"] = NA
datavect <- 1:20
#lm(y ~ x)
ci.lines.fun(x,y0,datavect,col.line="gray")
#
y1 <-co2nr.de1
y1[is.na(y1) | y1=="Inf"] = NA
points(x, y1, pch=16,  col =alpha("orange", 0.3))
ci.lines.fun(x,y1,datavect,col.line="orange")
#
y2 <-co2nr.de2
y2[is.na(y2) | y2=="Inf"] = NA
points(x, y2, pch=1,  col =alpha(1, 0.3))

ci.lines.fun(x,y2,datavect,col.line=1)
ci.lines.fun(x,y0,datavect,col.line="gray")
ci.lines.fun(x,y1,datavect,col.line="orange")


legend("topright", "e)", bty = "n", cex = 1.5)

#plot Co2 residuals vs SWC
plot(swc.hour, co2h.twe.resn,  ylim = c(-1,7),   
     xaxt="n",
     cex.axis=1.3, cex.lab=1.3, pch=15, col = alpha("gray", 0.3),
     xlab = expression("SWC"[10]~"(%)"),
     ylab = expression(norm.~CO[2]~residuals))
axis(1, at=seq(0,100,20),labels=T, cex.axis=1.3, cex.lab=1.3, col = 1)#, col.axis= "blue")
#add mean trend and confidence lines
x <- swc.hour
datavect <- 1:100
x[is.na(x) | x=="Inf"] = NA
y0 <-co2h.twe.resn
y0[is.na(y0) | y0=="Inf"] = NA

ci.lines.fun(x,y0,datavect,col.line="gray")
#
y1 <-co2nr.de1
y1[is.na(y1) | y1=="Inf"] = NA
points(x, y1, pch=16,  col =alpha("orange", 0.3))
ci.lines.fun(x,y1,datavect,col.line="orange")
#
y2 <-co2nr.de2
y2[is.na(y2) | y2=="Inf"] = NA
points(x, y2, pch=1,  col =alpha(1, 0.3))

ci.lines.fun(x,y2,datavect,col.line=1)
ci.lines.fun(x,y0,datavect,col.line="gray")
ci.lines.fun(x,y1,datavect,col.line="orange")
legend("topright", "f)", bty = "n", cex = 1.5)
legend("topleft", legend=c(expression(Y07.~xi~TW), 
                           expression(Y07.~xi~D~"p("~theta~"|SOC)"),
                           expression(Y07.~xi~D~"p("~theta~"|SOC-CO2)")),
       pch = c(15,16,1), lty= 1, lwd=2,col = c("gray","orange",1), bg = "white", bty = "n", cex=1.3, pt.cex = c(1.3,1.3,1.3))

