#
#Functions to calculate litter

# foliar litter (needles)
# Agren et al 2007 
# Pine : Turnover rate = 1.656 - 0.0231*Latitude (degres)
# Spruce : Turnover rate =  0.489 - 0.0063*Latitude (degres)
foliage.litter.swe <- function(Mf, spec, lat) {
  if (spec==1 ) {
    return(Mf*(1.656 - 0.0231*lat))} #1.656 - 0.0231*c(57:70)
  if (spec==2){
    return(Mf*(0.489 - 0.0063*lat))} #0.489 - 0.0063*c(57:70)
  if (spec==3) {
    return(Mf*0.79)} # same as for Finland, lowering biomas for the resorbed nutrients and carbohydrates during autumnal yellowing process
}


# Liski et al. 2006 Tab.1 biomas turnover rates (year-1) used to estimate the litter production of trees and ground vegetation
# branch litter
branch.litter <- function(Mb, spec) {
  if (spec==1) {
    return(Mb*0.02)}
  if (spec==2) {
    return(Mb*0.0125)}
  if (spec==3) {
    return(Mb*0.0135)}
}

# root litter
root.litter <- function(Mr, spec) {
  if (spec==1) {
    return(Mr*0.0184)}
  if (spec==2) {
    return(Mr*0.0125)}
  if (spec==3) {
    return(Mr*0.0135)}
}


# stump bark litter
stump.litter <- function(Mst, spec) {
  if (spec==1) {
    return(Mst*0.0029)}
  if (spec==2) {
    return(0)}
  if (spec==3) {
    return(Mst*0.0001)}
}

#reproductive origins  and stem bark
rep.stem.litter <- function(Ms, spec) {
  if (spec==1) {
    return(Ms*0.0052)}
  if (spec==2) {
    return(Ms*0.0027)}
  if (spec==3) {
    return(Ms*0.0029)}
}


#proportion of maximum temperature(1=max,0=min)
#temp.ratio <- (se.carbon$Tair.mean+abs(min(se.carbon$Tair.mean)))/
#                              diff(range(se.carbon$Tair.mean))
#se.carbon$temp.ratio<-temp.ratio
#plot(se.carbon$Tair.mean,temp.ratio)

# fine root litter with correction for temperature
fineroot.litter <- function(Mf,spec,temp.ratio) {#proportion of maximum temperature(1=max,0=min)
  #max differences for fine.root.turnover rate as related to south-north temperature difference
  dtr.s <-0.868-0.5
  dtr.p <-0.811-0.5
  dtr.b <-1-0.5
  
  if (spec==1) {
    return(Mf*(0.5+dtr.p*temp.ratio))} #uncorrected Mf*0.868
  if (spec==2) {
    return(Mf*(0.5+dtr.s*temp.ratio))} #uncorrected Mf*0.811
  if (spec==3) {
    return(Mf*(0.5+dtr.b*temp.ratio))} #uncorrected Mf*1
}




# Understorey litter ratios ##############################################


dwarfshrub.litter <- function(mass, comp) {
  if (comp=='abv'){
    return(0.37*mass) }
  if (comp=='bel'){
    return( 0.08767647*mass) }
  
}

# here calculating bel litter
# Varpujen hienojuurien osuus maanalaisesta kokonaisbiomassasta on keskimäärin 7 %.
# 1/1.7 = 0.5882353
# maavarret 20v 0.05
# tr (93*0.05+7*0.5882353) / 100 = 0.08767647

bryof.litter <- function(mass) {
  return(0.417*mass)
}

lichen.litter <- function(mass) {
  return(0.1*mass)
}

grass.litter <- function(mass,comp) {
  if (comp=='abv'){
    return(0.33*mass) }
  if (comp=='bel'){
    return((1/1.7)*mass) }
}

herb.litter <- function(mass,comp){
  if (comp=='abv'){
    return(1*mass) }
  if (comp=='bel'){
    return((1/1.7)*mass) }
}

