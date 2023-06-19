
##################################################################################################################
# Reconstruct figures of Tupek et al. 2023: 
# Modeling boreal forest’s mineral soil and peat C stock dynamics with Yasso07 model coupled with updated moisture modifier
#
# Run Yasso07 soil carbon model coupled with soil temp and moisture functions 
# calibrated with data from Vatiharju-Lakkasuo boreal frorest - mire ecotone in Finland 
# (Tupek et al. 2008, Tupek et al. 2015)
###################################################################################################################

# Boris Tupek
# Natural Resources Institute Finland (LUKE)
# June 2023
# boris.tupek@luke.fi


# CODE STRUCTURE:
# Figures 2,3
# 1) plot figures
#    based on the daily CO2 measurements, soil temperature and moisture, 
#    soil C stocks and litter data from 9 sites of forest-mire ecotone

# Figures 4 and 5 
# 2) run Yasso07 model in original 
# 3) run Yasso07 model with modified envi functions


rm(list=ls())

## preparation objects for plotting the figures#########
## color palettes
#install.packages('squash', dependencies = T)
library(squash)
#define colors
blackwhite.palette <- colorRampPalette(c("white", "black"), space = "rgb")
bw.colors <- blackwhite.palette(110)[10:110]
bw.palette <- colorRampPalette(c("black", "white"), space = "rgb")
bluegreen.palette <- colorRampPalette(c("blue", "green"), space = "rgb")

#A function to add arrows on the chart ###########
error.bar <- function(x, y, upper, lower, length=0.1,...){ #lower=upper, 
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
#############################################

# Rh,T,SWC data #################
r <- read.delim("D:/LUKE/PAPER-BG20/SUBMISSION/GMD23_YassoSWC/CODE-GMD/DATA/r.notveg_t5sm10_eco.g.ft_19.06.23.csv", header = T, skip = 0, sep=";")


## SOC interpolated for 1 m (FIG S1), and SOC of profiles for interpolation 
soilc.bg15 <- read.delim( "D:/LUKE/PAPER-BG20/SUBMISSION/GMD23_YassoSWC/CODE-GMD/DATA/soilc.bg15_19.06.23.csv", header = T, skip = 0, sep=";")
soilc.oct2015 <- read.delim( "D:/LUKE/PAPER-BG20/SUBMISSION/GMD23_YassoSWC/CODE-GMD/DATA/soilc.oct2015_19.06.23.csv", header = T, skip = 0, sep=";")
tSOC <- read.delim( "D:/LUKE/PAPER-BG20/SUBMISSION/GMD23_YassoSWC/CODE-GMD/DATA/tSOC_19.06.23.csv")$x
tSOCw <- read.delim("D:/LUKE/PAPER-BG20/SUBMISSION/GMD23_YassoSWC/CODE-GMD/DATA/tSOCw_19.06.23.csv")$x
tSOC.ci1 <- read.delim( "D:/LUKE/PAPER-BG20/SUBMISSION/GMD23_YassoSWC/CODE-GMD/DATA/tSOC.ci1_19.06.23.csv")$x
tSOC.ci99<- read.delim("D:/LUKE/PAPER-BG20/SUBMISSION/GMD23_YassoSWC/CODE-GMD/DATA/tSOC.ci99.06.23.csv")$x

ls()

#plot Fig 2 ##
source("D:/LUKE/PAPER-BG20/SUBMISSION/GMD23_YassoSWC/CODE-GMD/SCRIPTS/y7swc_gmd23_fig2.r")

# plot Fig 3 ##
source("D:/LUKE/PAPER-BG20/SUBMISSION/GMD23_YassoSWC/CODE-GMD/SCRIPTS/y7swc_gmd23_fig3.r")


##################################################################################################################
## SOC CO2 MODELLING #############################################################################################

#source yasso07 model function (correct version with litter input to all A,W,E,N pools compatible with SoilR, not just A which is erroneous)
source("D:/LUKE/PAPER-BG20/SCRIPTS-BG/yasso07Model.soilr.fi_function.r")
#this includes function for equilibrium
#list functions
lsf.str()

