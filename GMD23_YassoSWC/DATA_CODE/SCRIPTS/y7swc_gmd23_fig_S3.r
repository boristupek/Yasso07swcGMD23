
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
#LITTER INPUT  (plot Fig S2 and Fig S3)############
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

## Plot FIGS3 LITTER TIME SERIES ################################################################################################################
g9.clrs <- g.clrs(11)
ltot <- rowSums(litter.fc.month.awen_tot[, c("A.tot", "W.tot", "E.tot", "N.tot")])/10000 #convert from ha to m2
p.ix <- floor(seq(1,432/9, length.out=11))

par(mfrow=c(1,1), mar=c(5,2,1,1), oma=c(0,3,0,0))
for(i in 1:9){
  ix <- which(litter.fc.month.awen_tot$plot==i)
  
  date1 <- paste(litter.fc.month.awen_tot[ix,"year"],
                 litter.fc.month.awen_tot[ix,"month"],
                 11, sep ="-" )[p.ix]
  
  xlab.date1 <- format(as.POSIXct(as.Date(date1)),format="%b %y")
  
  if(i == 1){
    plot(ltot[ix], xaxt = "n",
         col = g9.clrs[(i+2)],
         ylab="litter",#expression(Litter~~group("(",kg~C~m^{-2}~month^{-1},")")),
         xlab = "time (months)",
         type ="l", log="y", ylim =range(ltot) )
    points(p.ix, ltot[ix][p.ix], pch = i,col = g9.clrs[(i+2)])
    axis(1, at=p.ix,labels=xlab.date1,cex.axis=1)
    
  }else{
    lines(ltot[ix], 
          col = g9.clrs[(i+2)])
    points(p.ix, ltot[ix][p.ix], pch = i,col = g9.clrs[(i+2)])
  }
}
abline(v=p.ix[5], col = 2, lty = 2)
legend("topright", c("CT","VT","MT","OMT","OMT+","KGK","KR","VSR1","VSR2"),
       ncol =3, col = g9.clrs[3:11], lty=1, pch = 1:9,
       bty="n", cex= 1)
text(p.ix[5],0.0001, "trenching", col = 2)
mtext(expression(Litter~~group("(",kg~C~m^{-2}~month^{-1},")")), 2, 0.7, outer=T, cex= 1.2)

######################################################################################################################