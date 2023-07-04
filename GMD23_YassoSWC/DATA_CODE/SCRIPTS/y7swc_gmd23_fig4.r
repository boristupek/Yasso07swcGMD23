
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

## Plot FIG 4 ##################################################################################################
## Environmental modifiers of soil temperature and moisture ##############################################

#separate scale
colors <- rev(greenblue.palette(110)[10:110])#heat.colors(901)

#CLIMATE SPACE for modifiers two directional changes (temperature, moisture) ######
t90=seq(-1,20,length.out=90)
w90=seq(1,90,length.out=90)/100
t90rep=rep(-1:20,each=90)
w90rep=rep(seq(0,90,length.out=90)/100,
           times=length(-1:20))

library(doBy)
names(r)
t5sm10m <- summaryBy(cbind(t5mo1, sm10mo1) ~ ft+month, data=r, FUN=mean) #function(x){mean(x,na.rm=T)}


#T, SWC functions separately and combined
temp.curve.fun <- function(q){
  q^((t90 -10)/10)
}

swc.curve.fun <- function(d,wopt){
  d^((wopt - w90*100)^2)
}

tw.climspace.fun <- function(q,d,wopt,wrep,trep){
  d^((wopt - wrep*100)^2)*q^((trep -10)/10)
  
}

d.soc<-map.soc["d"]
w.soc<-map.soc["wopt"]

d.socco2<-map.socco2["d"]
w.socco2<-map.socco2["wopt"]

q.soc<-as.numeric(map.soc["q"])
q.socco2<-map.socco2["q"]

map.soc
#d         q      wopt 
#0.9994776 4.3382941 7.0976308 

map.socco2
#d         q      wopt 
#0.9994907 3.4252899 5.0092200 

#climate space 1D LINES for 1 predicted modifier (W or T)
ew.soc <- swc.curve.fun(d.soc,w.soc)
ew.socco2 <- swc.curve.fun(d.socco2,w.socco2)

et.soc <- temp.curve.fun(q.soc)
et.socco2 <- temp.curve.fun(q.socco2)

#climate space 2D MATRICES of W and T predicted modifiers combined
e.soc <-tw.climspace.fun(q.soc,d.soc,w.soc,w90rep,t90rep)
me.soc = t(matrix(e.soc, nrow=90))

e.socco2 <-tw.climspace.fun(q.socco2,d.socco2,w.socco2,w90rep,t90rep)
me.socco2 = t(matrix(e.socco2, nrow=90)) 

#MS FIGURE 4 ## MODIFIERS #########

#figure
#modifiers in temp and swc climate space
#separate scale
colors <- rev(greenblue.palette(110)[10:110])#heat.colors(901)

par(mfrow=c(2,2), mar = c(5,5,3,2),oma=c(0,1,0,1))
#ft
x=seq(-1,20, length.out = 90)
plot(x,1:90, col="white", xlab= expression("T"[5]~"(°C)"),#"Soil Temperature (\u00B0C)",
     ylab=expression(xi[D]~"="~italic(f)~"(T"[5]~")"),
     main = "",type="l", ylim = c(0,5.1)) 
lines(x,et.soc, col = 1)
lines(x,et.socco2 , col = 2, lty=2)
legend("topleft",c(expression(xi[D]~~"~ p ("~ theta ~"| SOC )"), expression(xi[D]~~"~ p ("~ theta ~"| SOC-CO2 )")),
       lty=c(1, 2,1,1), col=c(1, 2,2,4), bty="n")
legend("bottomright","a)",bty = "n",cex=1.7)

#fswc
plot(1:90, col="white", xlab=expression("SWC"[10]~"(%)"),
     ylab=expression(xi[D]~"="~italic(f)~"(SWC"[10]~")"), 
     main = "",type="l", ylim = c(0,1.0)) 
lines(ew.soc, col = 1)
lines(ew.socco2, col = 2, lty=2)
legend("bottomright","b)",bty = "n",cex=1.7)

image(y=-1:20,x=seq(0,90,length.out=90), t(me.soc), col= colors,
      xlab="",#expression("Soil Temperature"[5]~"(°C)"),
      ylab="",#expression("Soil Water Content "[10]~~"(%)"),
      main = expression(xi[D]~"="~italic(f)~"(T"[5]~", SWC"[10]~")"~"~ p("~ theta ~"| SOC )"))
points(t5sm10m$sm10mo1.mean, t5sm10m$t5mo1.mean, pch=1, col ="white" ) #monthly corresponding temp5 and swc10
points( swc.9.mean, tn5.9.mean, pch=1, col =1)
contour(y=-1:20,x=seq(0,90,length.out=90),t(me.soc),add=TRUE,nlevels = 11)
contour(y=-1:20,x=seq(0,90,length.out=90),t(me.soc),add=TRUE,
        levels = c(0.05, 0.10, 0.20, 0.3))
contour(y=-1:20,x=seq(0,90,length.out=90),t(me.soc),add=TRUE,
        levels = 1, lwd=2)
mtext(expression("T"[5]~"(°C)"), side = 2, line = 2.3, outer = F, cex = 0.8)
mtext(expression("SWC"[10]~~"(%)"), side = 1, line = 2.7, outer = F, cex = 0.8)
legend("bottomright","c)",bty = "n",cex=1.7)

image(y=-1:20,x=seq(0,90,length.out=90), t(me.socco2), col= colors,
      xlab="",#expression("Soil Temperature"[5]~"(°C)"),
      ylab="",#expression("Soil Water Content "[10]~~"(%)"),
      main = expression(xi[D]~"="~italic(f)~"(T"[5]~", SWC"[10]~")"~"~ p("~theta ~"| SOC-CO2 )"))
points(t5sm10m$sm10mo1.mean, t5sm10m$t5mo1.mean, pch=1, col ="white" ) #monthly corresponding temp5 and swc10
points( swc.9.mean, tn5.9.mean, pch=1, col =1)
contour(y=-1:20,x=seq(0,90,length.out=90),t(me.socco2),add=TRUE,nlevels = 11)
contour(y=-1:20,x=seq(0,90,length.out=90),t(me.socco2),add=TRUE,
        levels = c(0.05, 0.10))
contour(y=-1:20,x=seq(0,90,length.out=90),t(me.socco2),add=TRUE,
        levels = 1, lwd=2)
mtext(expression("T"[5]~"(°C)"), side = 2, line = 2.3, outer = F, cex = 0.8)
mtext(expression("SWC"[10]~~"(%)"), side = 1, line = 2.7, outer = F, cex = 0.8)
legend("bottomright","d)",bty = "n",cex=1.7)
