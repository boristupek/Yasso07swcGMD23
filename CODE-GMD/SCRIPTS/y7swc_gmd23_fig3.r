
##################################################################################################################
# Reconstruct figures of Tupek et al. 2023: 
# Modeling boreal forest’s mineral soil and peat C stock dynamics with Yasso07 model coupled with updated moisture modifier
#
###################################################################################################################

# Boris Tupek
# Natural Resources Institute Finland (LUKE)
# June 2023
# boris.tupek@luke.fi

#note: standalone requires r data
# Rh,T,SWC data #################
#r <- read.delim("D:/LUKE/PAPER-BG20/SUBMISSION/GMD23_YassoSWC/CODE-GMD/DATA/r.notveg_t5sm10_eco.g.ft_19.06.23.csv", header = T, skip = 0, sep=";")

## FIG 3 # CO2 R +METEO seasonal ######################################################################################

#functions for fig 3
# extract copy/edit from the code below
plot.rn.9 <- function(i) {
  #for(i in 1:9){
  #i=1
  si <-which(r$ft ==i)
  ft.name <-c("CT","VT","MT","OMT","OMT+","KgK","KR","VSR1","VSR2")[i]
  if(i==1){
    plot(r$doy[si],r$rn[si], pch=16,ylab= "", 
         xaxt ='n',
         xlab="", #
         main="", #
         col= "white",
         xlim=c(90,350),
         ylim=c(-0.1,2.1) )
    axis(1, at=seq(0,365,50),labels=F) 
  } else{
    plot(r$doy[si],r$rn[si], pch=16,ylab= "", 
         xaxt ='n',
         xlab="", #
         main="", #
         col= "white",
         xlim=c(90,350),
         ylim=c(-0.1,2.1),
         yaxt="n")
    
    axis(1, at=seq(0,365,50),labels=F)
    ticks= seq(-0.1, 2.1,0.5)
    axis(2, at = ticks, labels = F)
  }
  
  for (j in 1:3){
    #j=1
    y=c(2004:2006)[j]
    sy = which(r$ft ==i & r$year==y)
    points(r$doy[sy],r$rn[sy],
           pch=c(16,1,16)[j],
           col = bw.palette(5)[j])
  }   
  mtext(outer=F, side=3, line=1, cex=1.1,text=ft.name)
  #}
}

plot.t5.9 <- function(i) {
  #par(mfrow=c(3,3))
  #for(i in 1:9){
  #i=1
  si <-which(r$ft ==i)
  ft.name <-c("CT","VT","MT","OMT","OMT+","KgK","KR","VSR1","VSR2")[i]
  if(i==1){
    plot(r$doy[si],r$t5mo1[si], pch=16,ylab= "", 
         xaxt ='n',
         xlab="", #
         main="", #
         col= "white",
         xlim=c(90,350),
         ylim=c(0,22) )
    
    axis(1, at=seq(0,365,50),labels=F) 
  } else{
    plot(r$doy[si],r$t5mo1[si], pch=16,ylab= "", 
         xaxt ='n',
         xlab="", #
         main="", #
         col= "white",
         xlim=c(90,350),
         ylim=c(0,22),
         yaxt="n")
    
    axis(1, at=seq(0,365,50),labels=F)
    ticks= seq(0, 22,5)
    axis(2, at = ticks, labels = F)
  }
  
  for (j in 1:3){
    #j=1
    y=c(2004:2006)[j]
    sy = which(r$ft ==i & r$year==y)
    points(r$doy[sy],r$t5mo1[sy], 
           pch=c(16,1,16)[j],
           col = bw.palette(5)[j])
    
  }   
  
}


plot.sm10mo1.9 <- function(i) {
  #for(i in 1:9){
  #i=1
  si <-which(r$ft ==i)
  ft.name <-c("CT","VT","MT","OMT","OMT+","KgK","KR","VSR1","VSR2")[i]
  if(i==1){
    plot(r$doy[si],r$sm10mo1[si], pch=16,ylab= "", 
         xaxt ='n',log="y",
         xlab="", #
         main="", #
         col= "white",
         xlim=c(90,350),
         ylim= c(5,101))
    
    axis(1, at=seq(0,365,50),labels=T) 
  } else{
    plot(r$doy[si],r$sm10mo1[si], pch=16,ylab= "", 
         xaxt ='n',log="y",
         xlab="", #
         main="", #
         col= "white",
         xlim=c(90,350),
         ylim=c(5,101),
         yaxt ='n')
    
    axis(1, at=seq(0,365,50),labels=T)
    ticks= c(5,10,20,50,100)
    axis(2, at = ticks, labels = F)
  }
  
  for (j in 1:3){
    #j=1
    y=c(2004:2006)[j]
    sy = which(r$ft ==i & r$year==y)
    points(r$doy[sy],r$sm10mo1[sy],
           pch=c(16,1,16)[j],
           col = bw.palette(5)[j])
  }
}

# PLOT fig 2 # CO2 R +METEO seasonal ############
# forest/mire types separately
#seasonal.Rco2.gm2h.T5.SWC10.forest.mire.types.H

par(mfcol=c(3,9))
par(mar=c(0,0,0,0), oma=c(6,6,3,1))

for (i in 1:9){
  #ft=1
  #respiration
  plot.rn.9(i)
  
  if(i == 1){
    legend(100,2.1,bty="n", title = "Year:",
           legend=c("wet","typical","dry" ),
           pch=c(16,1,16), col=bw.palette(5)[1:3], cex = 1.3)
    mtext(c("a)"), side = 3, line = -1,  at = 1,
          cex = 1.3)
    mtext(text=expression(R[h]~~group("(",gCO[2]~~m^2~~h^{-1},")")),
          side=2, line=3, cex=1.2)
  }
  
  #temperature
  plot.t5.9(i)
  if(i == 1){
    mtext(c("b)"), side = 3, line = -1,  at = 1,
          cex = 1.3)
    mtext(expression(T[5]~~"(°C)"),side=2, line = 3, cex=1.2)
  }
  #moisture
  plot.sm10mo1.9(i)
  if(i == 1){
    legend("topleft",bty="n", title = "Year:",
           legend=c("wet","typical","dry" ),
           pch=c(16,1,16), col=bw.palette(5)[1:3], cex = 1.3)
    mtext(c("c)"), side = 3, line = -1,  at = 1,
          cex = 1.3)
    mtext(expression(SWC[10]~~"(%)"),side=2, line = 3, cex=1.2)
  }
  
}

mtext(outer=T, side=1, line=3, cex=1.3,  text="Day of Year") #caption for x-axis

#################################################################################
