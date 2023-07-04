
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
#LITTER INPUT  (data to plot Fig S2 and Fig S3)############
#
#load modeled biomass and litter by method of finnish GHG described in Tupek et al 2019
#load(file=paste(paste0(path.data,"biomass.litter_trees.understory_years_12.02.20.RData"),sep=""))
#includes "ul" understory,"tl"- trees, "tl.msd" agregated for species, "tul.sum" agre. for sites,
#"tul.bp" for plotting (NOTE tree finerots excluded from the litte input as the tree roots were cut)

#load litter input from trenching for 2004 JULY Time series input!!!
#load(paste(paste0(path.data,"litter.trench.cut_coarse.fineroots_12.02.20.RData"),sep=""))
#includes "tree.w.coarseroot.cut04.awen.sum", "treeunder.nw.fineroot.cut04.awen"

#load lists of simulated AWEN data woody and nowoody litter 
#load(file=paste(paste0(path.data,"litter.eco_equilibrium_months_12.02.20.RData"),sep=""))
#includes lists "litter.eq.awen", "litter.natural.month.awen","litter.fc.month.awen"

## Plot FIGS2 ##################################################################################################
## Total LITTER input, STACKBAR components of tree and understory ##############################################

par(mfrow=c(3,1), mar=c(1,2,2,1), oma=c(5,5,0,0), xpd=TRUE)
x <- t(tul.bp.eq)/10000 #kg/ha convert to kg m2!!!!
#matrix in gray shade#tree litter
x.mids <- barplot(x, width =rep(0.5,9), col=g5.eq.clrs , axisnames=F, ylim = c(0,0.5), #axes=F, 
                  space=c(0.1,0,0,0, 0.2,0,0, 0.2,0), #col = b3.clrs,
                  las = 1, density = NULL, cex.axis = 1.3)
axis(1, at=seq(0.3,4.5,length.out = 9),labels= F, #c("CT","VT","MT","OMT","OMT+","KGK","KR","VSR1","VSR2"), 
     cex.axis=1.3)
x1 <- rbind(x[1,]/2,x[1,]/2) #understory total #add density lines
#understory half aboveground and half belowground
barplot(x1,col=1, axisnames=F, density=15, angle = c(120,30), 
        add=T, axes=F, width =rep(0.5,9),
        c(x.mids[1],x.mids[9],9), space=c(0.1,0,0,0, 0.2,0,0, 0.2,0))
mtext(expression("Litter (kg C m"^{-2}~"year"^{-1}~")"), 2, 2.9, cex = 1.2, outer =F)

xt <- tul.sum$totlit.eq.lu/10000 #convert to kg m2!!!!
xt.sd <- 2*tul.sum$totlit.rh.eq.sd.sum/10000+0.05*tul.sum$totlit.u/10000  #convert to kg m2!!!!

#add error bars (only for totals)
x_err.lower <- as.numeric(xt.sd)#/2
x_err.upper <- as.numeric(xt.sd)
error.bar(x.mids,xt, x_err.upper, x_err.lower)

legend(0.1,0.55,
       legend = c("foliar", "branch","stem+stump", "coarseroot","fineroot"), 
       fill = g5.eq.clrs[c(2:6,1,1)], 
       border= 1,
       bty = "n", ncol = 2, cex = 1.3, title = "Equilibrium:")

legend(1.5,0.55,
       legend = c("understory aboveground","understory belowground"),
       fill = T,
       col=g5.clrs[1],
       border= 1, 
       density=17,
       angle=c(30,120),
       bty = "n", ncol = 1, cex = 1.3, title = "")

## litter field campaign tree understory / LITTER PER MONTH!!!!!!!!!!!
#1) exclude fineroot and coarseroot litter, and understory litter
tul.bp.fc <- tul.bp.eq[,c("lit.needle.m.sum", "lit.branch.m.sum", "lit.stemstump.m.sum")]

x <- t(tul.bp.fc)/10000 #kg/ha convert to kg m2 MONTH!!!!

#total litter during field campaign remove cut litter from total litter
names(tul.bp.eq)
totlit.fc.lu <- tul.sum$totlit.eq.lu - rowSums(tul.bp.eq[,c("totlit.u", "lit.coarseroot.m.sum", "lit.fineroot.m.sum")])
#reduce sd
sd.fc.reduction <- totlit.fc.lu/tul.sum$totlit.eq.lu

#2) add dead root biomass from cutting it off by saw and preventing ingrowth by metalic collar

coarseroot.biom.cut04 <- rowSums(litter.trenchcut04[[1]][,3:6])
fineroot.biom.cut04 <- rowSums(litter.trenchcut04[[2]][,3:6])

x.cut <- rbind(t(tul.bp.fc)/12,coarseroot.biom.cut04,fineroot.biom.cut04)/10000 #kg/ha convert to kg m2!!!!

#tree cut roots biomass/litter one time input in July 2004
x.mids <- barplot(x.cut, width =rep(0.5,9),  #log="y",
                  col=g5.eq.clrs[2:6] , axisnames=F, ylim = c(0,2.5), #axes=F, 
                  space=c(0.1,0,0,0, 0.2,0,0, 0.2,0), #col = b3.clrs,
                  las = 1, density = NULL,
                  cex.axis = 1.2)
axis(1, at=seq(0.3,4.5,length.out = 9),labels=F,cex.axis=1.3)
mtext(expression("Litter (kg C m"^{-2}~"month"^{-1}~")"), 2, 2.9, cex = 1.2, outer =F)

legend(3.1,2.6,
       legend = c("foliar", "branch","stem+stump", "coarseroot","fineroot"), 
       fill = g5.eq.clrs[c(2:6,1,1)], 
       border= 1,
       bty = "n", ncol = 2, cex = 1.3, title = "July 2004 trenching, harvesting:")

#tree litter for time series
x.mids <- barplot(x, width =rep(0.5,9), col=g5.eq.clrs[2:4] , axisnames=F, ylim = c(0,0.15), #axes=F, 
                  space=c(0.1,0,0,0, 0.2,0,0, 0.2,0), #col = b3.clrs,
                  las = 1, density = NULL, cex.axis = 1.2)
axis(1, at=seq(0.3,4.5,length.out = 9),labels=c("CT","VT","MT","OMT","OMT+","KGK","KR","VSR1","VSR2"), 
     cex.axis=1.3)
mtext(expression("Litter (kg C m"^{-2}~"year"^{-1}~")"), 2, 3.3, cex = 1.2, outer =F)

xt <- totlit.fc.lu/10000 #convert to kg m2!!!!
xt.sd <- sd.fc.reduction*(2*tul.sum$totlit.rh.eq.sd.sum/10000+0.05*tul.sum$totlit.u/10000)  #convert to kg m2!!!!

#add error bars (only for totals)
x_err.lower <- as.numeric(xt.sd)#/2
x_err.upper <- as.numeric(xt.sd)
error.bar(x.mids,xt, x_err.upper, x_err.lower)

legend(3.1,0.15,
       legend = c("foliar", "branch","stem+stump"), 
       fill = g5.eq.clrs[2:4], 
       border= 1, 
       bty = "n", ncol = 1, cex = 1.3, title = "August 2004 onwards:")

mtext("Forest/mire types", 1, 3.3, cex = 1.3)
######################################################################################################################