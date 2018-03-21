##########################################################################################################################################
#                                                                                                                                        #
# Code to replicate the primary analysis of ecosystem multifunctionality as presented in Ramus et al. (2017) An invasive foundation      #
# species enhances multifunctionality in a coastal ecosystem. PNAS 114(32):8580-8585.                                                    #
#                                                                                                                                        #
# This code calculates multifunctionaity, fits candidate models using nonlinear least squares for each function individually, and        #
# generates the model selection table and corresponding figure presented in the paper.                                                   # 
#                                                                                                                                        #
# Developed by Aaron Ramus (aaron.ramus@gmail.com). Last updated 21 March 2018.                                                          #
#                                                                                                                                        #
##########################################################################################################################################


# Clear field
rm(list=ls(all.names=T))

# Load required libraries
library(AICcmodavg)
library(broom)
library(ggplot2)
library(multifunc)
library(nls2)
library(reshape2)

##########################################################################################################################################

# Read in data
meanPlots <- read.csv("1 mean plot-level responses.csv")

# Define ID variables
idVars <- c("Plot", "TrtPeg")

# Select predictor variable
predictor <- "Gcvr"
predictorName <- "Gracilaria cover (%)"

# Select response variables to be included in analysis
mfVars <- qw(
   Epi, 
   EpiRich, 
   #Dsln, 
   DslnFlip,
   Sed, 
   Nrsy, 
   NrsyRich, 
   Dcmp
   #Infa.sr, 
   #InfaRich.sr, 
   #Rays.sr, 
   #RaysFlip.sr,
   #Wfwl.sr
   ) 

# Calculate average multifunctionality (see Byrnes et al. 2014 Method Ecol Evol for usage of package 'multifunc')
meanPlots <- cbind(meanPlots, getStdAndMeanFunctions(meanPlots, mfVars))

# Create data frame for fits and plots
mfData <- meanPlots[,c(idVars, eval(predictor), mfVars, "meanFunction")]

# Calculate multifunctionality thresholds (again, see Byrnes et al. 2014 Method Ecol Evol for usage of package 'multifunc')
mfThresh <- getFuncsMaxed(meanPlots, mfVars, threshmin=0.10, threshmax=0.90, threshstep=.1, maxN=9, prepend=c(idVars, eval(predictor)))

# Add thresholds to mfData
mfThresh$thresholds <- as.integer(mfThresh$thresholds*100)
mfData$funcMaxed10 <- subset(mfThresh, thresholds==10)$funcMaxed
mfData$funcMaxed20 <- subset(mfThresh, thresholds==20)$funcMaxed
mfData$funcMaxed30 <- subset(mfThresh, thresholds==30)$funcMaxed
mfData$funcMaxed40 <- subset(mfThresh, thresholds==40)$funcMaxed
mfData$funcMaxed50 <- subset(mfThresh, thresholds==50)$funcMaxed
mfData$funcMaxed60 <- subset(mfThresh, thresholds==60)$funcMaxed
mfData$funcMaxed70 <- subset(mfThresh, thresholds==70)$funcMaxed
mfData$funcMaxed80 <- subset(mfThresh, thresholds==80)$funcMaxed
mfData$funcMaxed90 <- subset(mfThresh, thresholds==90)$funcMaxed

# Replace DslnFlip with unreflected dissolution for fits and plotting
mfData$DslnFlip <- meanPlots$Dsln
colnames(mfData)[6] <- "Dsln"

##########################################################################################################################################

# Create function to return information about the fitted nls model
getModFit <- function(x) {
   df.modfit <- data.frame(
      Model=as.character(paste(
         summary(x)$'formula'[2], 
         summary(x)$'formula'[1], 
         summary(x)$'formula'[3])), 
      a=round(tidy(x)$estimate[1], 3), 
      pa=tidy(x)$'p.value'[1], 
      b=tidy(x)$estimate[2], 
      pb=tidy(x)$'p.value'[2], 
      c=tidy(x)$estimate[3], 
      pc=tidy(x)$'p.value'[3], 
      RSE=glance(x)[1], 
      df=glance(x)[8], 
      isConv=glance(x)[2], 
      mTol=glance(x)[3], 
      LL2=glance(x)[4], 
      AIC=glance(x)[5], 
      BIC=glance(x)[6], 
      RSS=glance(x)[7], 
      TSS=sum((xy$y-mean(xy$y))^2),
      fvalue=anova(Ho, x)[2,5],
      pvalue=anova(Ho, x)[2,6]
   )
}

# Create parameter dataframes to fit models with package nls2
a <- data.frame(a=c(-100, -50, -10, -5, -1, -0.5, -0.001, 0, 0.001, 0.5, 1, 5, 10, 50, 100))
b <- data.frame(b=c(rep(-100, 15), rep(-50, 15),  rep(-10, 15),  rep(-5, 15),  rep(-1, 15), rep(-0.5, 15), rep(-0.001, 15),  rep(0, 15), rep(0.001, 15),  rep(0.5, 15), rep(1, 15), rep(5, 15), rep(10, 15), rep(50, 15), rep(100, 15)))      
ab <- data.frame(a=rep(a$a, 15), b=b$b)
c <- data.frame(c=c(rep(-100, 225), rep(-50, 225),  rep(-10, 225),  rep(-5, 225),  rep(-1, 225), rep(-0.5, 225), rep(-0.001, 225),  rep(0, 225), rep(0.001, 225),  rep(0.5, 225), rep(1, 225), rep(5, 225), rep(10, 225), rep(50, 225), rep(100, 225)))      
abc <- data.frame(a=rep(ab$a, 15), b=rep(ab$b, 15), c=c$c)
bc <- ab
ac <- ab
colnames(bc) <- c("b", "c")
colnames(ac) <- c("a", "c")

# Define candidate model formulas
HoFormula <- y ~ a  #null
LmFormula <- y ~ a + b * x  #linear
LogFormula <- y ~ a + b * log(x + 1)  #log
HypFormula <- y ~ a * x/(b + x)  #hyperbolic
PowFormula <- y ~ a + b * x^c  #power

# Create an object to stuff model selection tables into later
cum.nlsTab <- NULL

##########################################################################################################################################

# Fit each candidate model to each function individually. This first one is annotated but the rest are not.

# Select response variable to fit models
response <- "Epi"
responseName <- "Epifauna abundance (# m^-2)"

# Create XY data frame to fit models
xy <- mfData[,c(paste(predictor), paste(response))]
colnames(xy) <- c("x", "y")

# Quick plot
p <- qplot(x, y, data=xy) + labs(x=paste(eval(predictorName)), y=paste(eval(responseName)))
p

# Fit null model
HoStart <- nls2(HoFormula, data=xy, start=a, algorithm="brute-force")
Ho <- nls2(HoFormula, data=xy, start=HoStart)

# Plot fitted null model
pHo <- stat_function(fun=function(x){coef(Ho)[1]+0*x}, aes(col="Null"))
p+pHo

# Fit linear model
LmStart <- nls2(LmFormula, data=xy, start=ab, algorithm="brute-force")
Lm <- nls2(LmFormula, data=xy, start=LmStart)
pLm <-stat_function(fun=function(x){coef(Lm)[1]+coef(Lm)[2]*x}, aes(col="Linear"))
p+pLm

# Fit log model
LogStart <- nls2(LogFormula, data=xy, start=ab, algorithm="brute-force")
Log <- nls2(LogFormula, data=xy, start=LogStart)
pLog <- stat_function(fun=function(x){coef(Log)[1]+coef(Log)[2]*log(x+1)}, aes(col="Log"))
p+pLog

#  Fit hyperbolic model
HypStart <- nls2(HypFormula, data=xy, start=ab, algorithm="brute-force")
#Hyp <- nls2(HypFormula, data=xy, start=HypStart)
Hyp <- nls2(HypFormula, data=xy, start=c(a=75, b=1), nls.control(maxiter=200))  # can try to fit manually
pHyp <- stat_function(fun=function(x){coef(Hyp)[1]*x/(coef(Hyp)[2]+x)}, aes(col="Hyperbolic"))
p+pHyp

#  Fit power model
PowStart <- nls2(PowFormula, data=xy, start=abc, algorithm="brute-force")
Pow <- nls2(PowFormula, data=xy, start=PowStart)
#Pow <- nls2(PowFormula, data=xy, start=c(a=0.2, b=1, c=1))  # can try to fit manually
pPow <- stat_function(fun=function(x){coef(Pow)[1]+coef(Pow)[2]*x^coef(Pow)[3]}, aes(col="Power"))
p+pPow

# All the candidate models need to be fit succussfully before the below file is sourced 
# Otherwise you will need to create a second version of this file and modify it accordingly
source("3 nlsTab.R")

##########################################################################################################################################

response <- "EpiRich"
responseName <- "Epifauna richness (taxa m^-2)"

xy <- mfData[,c(paste(predictor), paste(response))]
colnames(xy) <- c("x", "y")
p <- qplot(x, y, data=xy) + labs(x=paste(eval(predictorName)), y=paste(eval(responseName)))
p

HoStart <- nls2(HoFormula, data=xy, start=a, algorithm="brute-force")
Ho <- nls2(HoFormula, data=xy, start=HoStart)
pHo <- stat_function(fun=function(x){coef(Ho)[1]+0*x}, aes(col="Null"))
p+pHo

LmStart <- nls2(LmFormula, data=xy, start=ab, algorithm="brute-force")
Lm <- nls2(LmFormula, data=xy, start=LmStart)
pLm <-stat_function(fun=function(x){coef(Lm)[1]+coef(Lm)[2]*x}, aes(col="Linear"))
p+pLm

LogStart <- nls2(LogFormula, data=xy, start=ab, algorithm="brute-force")
Log <- nls2(LogFormula, data=xy, start=LogStart)
pLog <- stat_function(fun=function(x){coef(Log)[1]+coef(Log)[2]*log(x+1)}, aes(col="Log"))
p+pLog

HypStart <- nls2(HypFormula, data=xy, start=ab, algorithm="brute-force")
Hyp <- nls2(HypFormula, data=xy, start=c(a=75, b=1), nls.control(maxiter=200))
pHyp <- stat_function(fun=function(x){coef(Hyp)[1]*x/(coef(Hyp)[2]+x)}, aes(col="Hyperbolic"))
p+pHyp

PowStart <- nls2(PowFormula, data=xy, start=abc, algorithm="brute-force")
Pow <- nls2(PowFormula, data=xy, start=PowStart)
pPow <- stat_function(fun=function(x){coef(Pow)[1]+coef(Pow)[2]*x^coef(Pow)[3]}, aes(col="Power"))
p+pPow

source("3 nlsTab.R")

##########################################################################################################################################

response <- "Dsln"
responseName <- "Dissolution (mass lost in g d^-1)"

xy <- mfData[,c(paste(predictor), paste(response))]
colnames(xy) <- c("x", "y")
p <- qplot(x, y, data=xy) + labs(x=paste(eval(predictorName)), y=paste(eval(responseName)))
p

HoStart <- nls2(HoFormula, data=xy, start=a, algorithm="brute-force")
Ho <- nls2(HoFormula, data=xy, start=HoStart)
pHo <- stat_function(fun=function(x){coef(Ho)[1]+0*x}, aes(col="Null"))
p+pHo

LmStart <- nls2(LmFormula, data=xy, start=ab, algorithm="brute-force")
Lm <- nls2(LmFormula, data=xy, start=LmStart)
pLm <-stat_function(fun=function(x){coef(Lm)[1]+coef(Lm)[2]*x}, aes(col="Linear"))
p+pLm

LogStart <- nls2(LogFormula, data=xy, start=ab, algorithm="brute-force")
Log <- nls2(LogFormula, data=xy, start=LogStart)
pLog <- stat_function(fun=function(x){coef(Log)[1]+coef(Log)[2]*log(x+1)}, aes(col="Log"))
p+pLog

HypStart <- nls2(HypFormula, data=xy, start=ab, algorithm="brute-force")
Hyp <- nls2(HypFormula, data=xy, start=HypStart)
pHyp <- stat_function(fun=function(x){coef(Hyp)[1]*x/(coef(Hyp)[2]+x)}, aes(col="Hyperbolic"))
p+pHyp

PowStart <- nls2(PowFormula, data=xy, start=abc, algorithm="brute-force")
Pow <- nls2(PowFormula, data=xy, start=c(a=40, b=1, c=1))
pPow <- stat_function(fun=function(x){coef(Pow)[1]+coef(Pow)[2]*x^coef(Pow)[3]}, aes(col="Power"))
p+pPow

source("3 nlsTab.R")

##########################################################################################################################################

response <- "Sed"
responseName <- "Sediment stabilization (delta cm mo^-1)"

xy <- mfData[,c(paste(predictor), paste(response))]
colnames(xy) <- c("x", "y")
p <- qplot(x, y, data=xy) + labs(x=paste(eval(predictorName)), y=paste(eval(responseName)))
p

HoStart <- nls2(HoFormula, data=xy, start=a, algorithm="brute-force")
Ho <- nls2(HoFormula, data=xy, start=HoStart)
pHo <- stat_function(fun=function(x){coef(Ho)[1]+0*x}, aes(col="Null"))
p+pHo

LmStart <- nls2(LmFormula, data=xy, start=ab, algorithm="brute-force")
Lm <- nls2(LmFormula, data=xy, start=LmStart)
pLm <-stat_function(fun=function(x){coef(Lm)[1]+coef(Lm)[2]*x}, aes(col="Linear"))
p+pLm

LogStart <- nls2(LogFormula, data=xy, start=ab, algorithm="brute-force")
Log <- nls2(LogFormula, data=xy, start=LogStart)
pLog <- stat_function(fun=function(x){coef(Log)[1]+coef(Log)[2]*log(x+1)}, aes(col="Log"))
p+pLog

HypStart <- nls2(HypFormula, data=xy, start=ab, algorithm="brute-force")
Hyp <- nls2(HypFormula, data=xy, start=c(a=0.1, b=-5), nls.control(maxiter=200))
pHyp <- stat_function(fun=function(x){coef(Hyp)[1]*x/(coef(Hyp)[2]+x)}, aes(col="Hyperbolic"))
p+pHyp

PowStart <- nls2(PowFormula, data=xy, start=abc, algorithm="brute-force")
Pow <- nls2(PowFormula, data=xy, start=c(a=1, b=1, c=2.5)) 
pPow <- stat_function(fun=function(x){coef(Pow)[1]+coef(Pow)[2]*x^coef(Pow)[3]}, aes(col="Power"))
p+pPow

source("3 nlsTab.R")

##########################################################################################################################################

response <- "Nrsy"
responseName <- "Nursery abundance (# m^-2)"

xy <- mfData[,c(paste(predictor), paste(response))]
colnames(xy) <- c("x", "y")
p <- qplot(x, y, data=xy) + labs(x=paste(eval(predictorName)), y=paste(eval(responseName)))
p

HoStart <- nls2(HoFormula, data=xy, start=a, algorithm="brute-force")
Ho <- nls2(HoFormula, data=xy, start=HoStart)
pHo <- stat_function(fun=function(x){coef(Ho)[1]+0*x}, aes(col="Null"))
p+pHo

LmStart <- nls2(LmFormula, data=xy, start=ab, algorithm="brute-force")
Lm <- nls2(LmFormula, data=xy, start=LmStart)
pLm <-stat_function(fun=function(x){coef(Lm)[1]+coef(Lm)[2]*x}, aes(col="Linear"))
p+pLm

LogStart <- nls2(LogFormula, data=xy, start=ab, algorithm="brute-force")
Log <- nls2(LogFormula, data=xy, start=LogStart)
pLog <- stat_function(fun=function(x){coef(Log)[1]+coef(Log)[2]*log(x+1)}, aes(col="Log"))
p+pLog

HypStart <- nls2(HypFormula, data=xy, start=ab, algorithm="brute-force")
Hyp <- nls2(HypFormula, data=xy, start=HypStart)
pHyp <- stat_function(fun=function(x){coef(Hyp)[1]*x/(coef(Hyp)[2]+x)}, aes(col="Hyperbolic"))
p+pHyp

PowStart <- nls2(PowFormula, data=xy, start=abc, algorithm="brute-force")
Pow <- nls2(PowFormula, data=xy, start=PowStart)
pPow <- stat_function(fun=function(x){coef(Pow)[1]+coef(Pow)[2]*x^coef(Pow)[3]}, aes(col="Power"))
p+pPow

source("3 nlsTab.R")

##########################################################################################################################################

response <- "NrsyRich"
responseName <- "Nursery richness (taxa m^-2)"

xy <- mfData[,c(paste(predictor), paste(response))]
colnames(xy) <- c("x", "y")
p <- qplot(x, y, data=xy) + labs(x=paste(eval(predictorName)), y=paste(eval(responseName)))
p

HoStart <- nls2(HoFormula, data=xy, start=a, algorithm="brute-force")
Ho <- nls2(HoFormula, data=xy, start=HoStart)
pHo <- stat_function(fun=function(x){coef(Ho)[1]+0*x}, aes(col="Null"))
p+pHo

LmStart <- nls2(LmFormula, data=xy, start=ab, algorithm="brute-force")
Lm <- nls2(LmFormula, data=xy, start=LmStart)
pLm <-stat_function(fun=function(x){coef(Lm)[1]+coef(Lm)[2]*x}, aes(col="Linear"))
p+pLm

LogStart <- nls2(LogFormula, data=xy, start=ab, algorithm="brute-force")
Log <- nls2(LogFormula, data=xy, start=LogStart)
pLog <- stat_function(fun=function(x){coef(Log)[1]+coef(Log)[2]*log(x+1)}, aes(col="Log"))
p+pLog

HypStart <- nls2(HypFormula, data=xy, start=ab, algorithm="brute-force")
Hyp <- nls2(HypFormula, data=xy, start=HypStart)
pHyp <- stat_function(fun=function(x){coef(Hyp)[1]*x/(coef(Hyp)[2]+x)}, aes(col="Hyperbolic"))
p+pHyp

PowStart <- nls2(PowFormula, data=xy, start=abc, algorithm="brute-force")
Pow <- nls2(PowFormula, data=xy, start=PowStart)
pPow <- stat_function(fun=function(x){coef(Pow)[1]+coef(Pow)[2]*x^coef(Pow)[3]}, aes(col="Power"))
p+pPow

source("3 nlsTab.R")

##########################################################################################################################################

response <- "Dcmp"
responseName <- "Decomposition (mass lost in g mo^-1)"

xy <- mfData[,c(paste(predictor), paste(response))]
colnames(xy) <- c("x", "y")
p <- qplot(x, y, data=xy) + labs(x=paste(eval(predictorName)), y=paste(eval(responseName)))
p

HoStart <- nls2(HoFormula, data=xy, start=a, algorithm="brute-force")
Ho <- nls2(HoFormula, data=xy, start=HoStart)
pHo <- stat_function(fun=function(x){coef(Ho)[1]+0*x}, aes(col="Null"))
p+pHo

LmStart <- nls2(LmFormula, data=xy, start=ab, algorithm="brute-force")
Lm <- nls2(LmFormula, data=xy, start=LmStart)
pLm <-stat_function(fun=function(x){coef(Lm)[1]+coef(Lm)[2]*x}, aes(col="Linear"))
p+pLm

LogStart <- nls2(LogFormula, data=xy, start=ab, algorithm="brute-force")
Log <- nls2(LogFormula, data=xy, start=LogStart)
pLog <- stat_function(fun=function(x){coef(Log)[1]+coef(Log)[2]*log(x+1)}, aes(col="Log"))
p+pLog

HypStart <- nls2(HypFormula, data=xy, start=ab, algorithm="brute-force")
Hyp <- nls2(HypFormula, data=xy, start=HypStart)
pHyp <- stat_function(fun=function(x){coef(Hyp)[1]*x/(coef(Hyp)[2]+x)}, aes(col="Hyperbolic"))
p+pHyp

PowStart <- nls2(PowFormula, data=xy, start=abc, algorithm="brute-force")
Pow <- nls2(PowFormula, data=xy, start=c(a=1, b=1, c=2.5)) 
pPow <- stat_function(fun=function(x){coef(Pow)[1]+coef(Pow)[2]*x^coef(Pow)[3]}, aes(col="Power"))
p+pPow

source("3 nlsTab.R")

##########################################################################################################################################

response <- "meanFunction"
responseName <- "Multifunctionality (%)"

xy <- mfData[,c(paste(predictor), paste(response))]
colnames(xy) <- c("x", "y")
p <- qplot(x, y, data=xy) + labs(x=paste(eval(predictorName)), y=paste(eval(responseName)))
p

HoStart <- nls2(HoFormula, data=xy, start=a, algorithm="brute-force")
Ho <- nls2(HoFormula, data=xy, start=HoStart)
pHo <- stat_function(fun=function(x){coef(Ho)[1]+0*x}, aes(col="Null"))
p+pHo

LmStart <- nls2(LmFormula, data=xy, start=ab, algorithm="brute-force")
Lm <- nls2(LmFormula, data=xy, start=LmStart)
pLm <-stat_function(fun=function(x){coef(Lm)[1]+coef(Lm)[2]*x}, aes(col="Linear"))
p+pLm

LogStart <- nls2(LogFormula, data=xy, start=ab, algorithm="brute-force")
Log <- nls2(LogFormula, data=xy, start=LogStart)
pLog <- stat_function(fun=function(x){coef(Log)[1]+coef(Log)[2]*log(x+1)}, aes(col="Log"))
p+pLog

HypStart <- nls2(HypFormula, data=xy, start=ab, algorithm="brute-force")
Hyp <- nls2(HypFormula, data=xy, start=HypStart)
pHyp <- stat_function(fun=function(x){coef(Hyp)[1]*x/(coef(Hyp)[2]+x)}, aes(col="Hyperbolic"))
p+pHyp

PowStart <- nls2(PowFormula, data=xy, start=abc, algorithm="brute-force")
Pow <- nls2(PowFormula, data=xy, start=c(a=0.2, b=1, c=1)) 
pPow <- stat_function(fun=function(x){coef(Pow)[1]+coef(Pow)[2]*x^coef(Pow)[3]}, aes(col="Power"))
p+pPow

source("3 nlsTab.R")

##########################################################################################################################################

response <- "funcMaxed10"
responseName <- "Number of functions above 10% threshold"

xy <- mfData[,c(paste(predictor), paste(response))]
colnames(xy) <- c("x", "y")
p <- qplot(x, y, data=xy) + labs(x=paste(eval(predictorName)), y=paste(eval(responseName)))
p

HoStart <- nls2(HoFormula, data=xy, start=a, algorithm="brute-force")
Ho <- nls2(HoFormula, data=xy, start=HoStart)
pHo <- stat_function(fun=function(x){coef(Ho)[1]+0*x}, aes(col="Null"))
p+pHo

LmStart <- nls2(LmFormula, data=xy, start=ab, algorithm="brute-force")
Lm <- nls2(LmFormula, data=xy, start=LmStart)
pLm <-stat_function(fun=function(x){coef(Lm)[1]+coef(Lm)[2]*x}, aes(col="Linear"))
p+pLm

LogStart <- nls2(LogFormula, data=xy, start=ab, algorithm="brute-force")
Log <- nls2(LogFormula, data=xy, start=LogStart)
pLog <- stat_function(fun=function(x){coef(Log)[1]+coef(Log)[2]*log(x+1)}, aes(col="Log"))
p+pLog

HypStart <- nls2(HypFormula, data=xy, start=ab, algorithm="brute-force")
Hyp <- nls2(HypFormula, data=xy, start=HypStart)
#Hyp <- nls2(HypFormula, data=xy, start=c(a=0.1, b=-5), nls.control(maxiter=200))
pHyp <- stat_function(fun=function(x){coef(Hyp)[1]*x/(coef(Hyp)[2]+x)}, aes(col="Hyperbolic"))
p+pHyp

PowStart <- nls2(PowFormula, data=xy, start=abc, algorithm="brute-force")
Pow <- nls2(PowFormula, data=xy, start=PowStart)
pPow <- stat_function(fun=function(x){coef(Pow)[1]+coef(Pow)[2]*x^coef(Pow)[3]}, aes(col="Power"))
p+pPow

source("3 nlsTab.R")

##########################################################################################################################################

response <- "funcMaxed20"
responseName <- "Number of functions above 20% threshold"

xy <- mfData[,c(paste(predictor), paste(response))]
colnames(xy) <- c("x", "y")
p <- qplot(x, y, data=xy) + labs(x=paste(eval(predictorName)), y=paste(eval(responseName)))
p

HoStart <- nls2(HoFormula, data=xy, start=a, algorithm="brute-force")
Ho <- nls2(HoFormula, data=xy, start=HoStart)
pHo <- stat_function(fun=function(x){coef(Ho)[1]+0*x}, aes(col="Null"))
p+pHo

LmStart <- nls2(LmFormula, data=xy, start=ab, algorithm="brute-force")
Lm <- nls2(LmFormula, data=xy, start=LmStart)
pLm <-stat_function(fun=function(x){coef(Lm)[1]+coef(Lm)[2]*x}, aes(col="Linear"))
p+pLm

LogStart <- nls2(LogFormula, data=xy, start=ab, algorithm="brute-force")
Log <- nls2(LogFormula, data=xy, start=LogStart)
pLog <- stat_function(fun=function(x){coef(Log)[1]+coef(Log)[2]*log(x+1)}, aes(col="Log"))
p+pLog

HypStart <- nls2(HypFormula, data=xy, start=ab, algorithm="brute-force")
Hyp <- nls2(HypFormula, data=xy, start=HypStart)
pHyp <- stat_function(fun=function(x){coef(Hyp)[1]*x/(coef(Hyp)[2]+x)}, aes(col="Hyperbolic"))
p+pHyp

PowStart <- nls2(PowFormula, data=xy, start=abc, algorithm="brute-force")
Pow <- nls2(PowFormula, data=xy, start=PowStart)
pPow <- stat_function(fun=function(x){coef(Pow)[1]+coef(Pow)[2]*x^coef(Pow)[3]}, aes(col="Power"))
p+pPow

source("3 nlsTab.R")

##########################################################################################################################################

response <- "funcMaxed30"
responseName <- "Number of functions above 30% threshold"

xy <- mfData[,c(paste(predictor), paste(response))]
colnames(xy) <- c("x", "y")
p <- qplot(x, y, data=xy) + labs(x=paste(eval(predictorName)), y=paste(eval(responseName)))
p

HoStart <- nls2(HoFormula, data=xy, start=a, algorithm="brute-force")
Ho <- nls2(HoFormula, data=xy, start=HoStart)
pHo <- stat_function(fun=function(x){coef(Ho)[1]+0*x}, aes(col="Null"))
p+pHo

LmStart <- nls2(LmFormula, data=xy, start=ab, algorithm="brute-force")
Lm <- nls2(LmFormula, data=xy, start=LmStart)
pLm <-stat_function(fun=function(x){coef(Lm)[1]+coef(Lm)[2]*x}, aes(col="Linear"))
p+pLm

LogStart <- nls2(LogFormula, data=xy, start=ab, algorithm="brute-force")
Log <- nls2(LogFormula, data=xy, start=LogStart)
pLog <- stat_function(fun=function(x){coef(Log)[1]+coef(Log)[2]*log(x+1)}, aes(col="Log"))
p+pLog

HypStart <- nls2(HypFormula, data=xy, start=ab, algorithm="brute-force")
Hyp <- nls2(HypFormula, data=xy, start=HypStart)
pHyp <- stat_function(fun=function(x){coef(Hyp)[1]*x/(coef(Hyp)[2]+x)}, aes(col="Hyperbolic"))
p+pHyp

PowStart <- nls2(PowFormula, data=xy, start=abc, algorithm="brute-force")
Pow <- nls2(PowFormula, data=xy, start=PowStart)
pPow <- stat_function(fun=function(x){coef(Pow)[1]+coef(Pow)[2]*x^coef(Pow)[3]}, aes(col="Power"))
p+pPow

source("3 nlsTab.R")

##########################################################################################################################################

response <- "funcMaxed40"
responseName <- "Number of functions above 40% threshold"

xy <- mfData[,c(paste(predictor), paste(response))]
colnames(xy) <- c("x", "y")
p <- qplot(x, y, data=xy) + labs(x=paste(eval(predictorName)), y=paste(eval(responseName)))
p

HoStart <- nls2(HoFormula, data=xy, start=a, algorithm="brute-force")
Ho <- nls2(HoFormula, data=xy, start=HoStart)
pHo <- stat_function(fun=function(x){coef(Ho)[1]+0*x}, aes(col="Null"))
p+pHo

LmStart <- nls2(LmFormula, data=xy, start=ab, algorithm="brute-force")
Lm <- nls2(LmFormula, data=xy, start=LmStart)
pLm <-stat_function(fun=function(x){coef(Lm)[1]+coef(Lm)[2]*x}, aes(col="Linear"))
p+pLm

LogStart <- nls2(LogFormula, data=xy, start=ab, algorithm="brute-force")
Log <- nls2(LogFormula, data=xy, start=LogStart)
pLog <- stat_function(fun=function(x){coef(Log)[1]+coef(Log)[2]*log(x+1)}, aes(col="Log"))
p+pLog

HypStart <- nls2(HypFormula, data=xy, start=ab, algorithm="brute-force")
Hyp <- nls2(HypFormula, data=xy, start=HypStart)
pHyp <- stat_function(fun=function(x){coef(Hyp)[1]*x/(coef(Hyp)[2]+x)}, aes(col="Hyperbolic"))
p+pHyp

PowStart <- nls2(PowFormula, data=xy, start=abc, algorithm="brute-force")
Pow <- nls2(PowFormula, data=xy, start=PowStart)
pPow <- stat_function(fun=function(x){coef(Pow)[1]+coef(Pow)[2]*x^coef(Pow)[3]}, aes(col="Power"))
p+pPow

source("3 nlsTab.R")

##########################################################################################################################################

response <- "funcMaxed50"
responseName <- "Number of functions above 50% threshold"

xy <- mfData[,c(paste(predictor), paste(response))]
colnames(xy) <- c("x", "y")
p <- qplot(x, y, data=xy) + labs(x=paste(eval(predictorName)), y=paste(eval(responseName)))
p

HoStart <- nls2(HoFormula, data=xy, start=a, algorithm="brute-force")
Ho <- nls2(HoFormula, data=xy, start=HoStart)
pHo <- stat_function(fun=function(x){coef(Ho)[1]+0*x}, aes(col="Null"))
p+pHo

LmStart <- nls2(LmFormula, data=xy, start=ab, algorithm="brute-force")
Lm <- nls2(LmFormula, data=xy, start=LmStart)
pLm <-stat_function(fun=function(x){coef(Lm)[1]+coef(Lm)[2]*x}, aes(col="Linear"))
p+pLm

LogStart <- nls2(LogFormula, data=xy, start=ab, algorithm="brute-force")
Log <- nls2(LogFormula, data=xy, start=LogStart)
pLog <- stat_function(fun=function(x){coef(Log)[1]+coef(Log)[2]*log(x+1)}, aes(col="Log"))
p+pLog

HypStart <- nls2(HypFormula, data=xy, start=ab, algorithm="brute-force")
Hyp <- nls2(HypFormula, data=xy, start=HypStart)
pHyp <- stat_function(fun=function(x){coef(Hyp)[1]*x/(coef(Hyp)[2]+x)}, aes(col="Hyperbolic"))
p+pHyp

PowStart <- nls2(PowFormula, data=xy, start=abc, algorithm="brute-force")
Pow <- nls2(PowFormula, data=xy, start=PowStart)
pPow <- stat_function(fun=function(x){coef(Pow)[1]+coef(Pow)[2]*x^coef(Pow)[3]}, aes(col="Power"))
p+pPow

source("3 nlsTab.R")

##########################################################################################################################################

response <- "funcMaxed60"
responseName <- "Number of functions above 60% threshold"

xy <- mfData[,c(paste(predictor), paste(response))]
colnames(xy) <- c("x", "y")
p <- qplot(x, y, data=xy) + labs(x=paste(eval(predictorName)), y=paste(eval(responseName)))
p

HoStart <- nls2(HoFormula, data=xy, start=a, algorithm="brute-force")
Ho <- nls2(HoFormula, data=xy, start=HoStart)
pHo <- stat_function(fun=function(x){coef(Ho)[1]+0*x}, aes(col="Null"))
p+pHo

LmStart <- nls2(LmFormula, data=xy, start=ab, algorithm="brute-force")
Lm <- nls2(LmFormula, data=xy, start=LmStart)
pLm <-stat_function(fun=function(x){coef(Lm)[1]+coef(Lm)[2]*x}, aes(col="Linear"))
p+pLm

LogStart <- nls2(LogFormula, data=xy, start=ab, algorithm="brute-force")
Log <- nls2(LogFormula, data=xy, start=LogStart)
pLog <- stat_function(fun=function(x){coef(Log)[1]+coef(Log)[2]*log(x+1)}, aes(col="Log"))
p+pLog

HypStart <- nls2(HypFormula, data=xy, start=ab, algorithm="brute-force")
Hyp <- nls2(HypFormula, data=xy, start=HypStart)
pHyp <- stat_function(fun=function(x){coef(Hyp)[1]*x/(coef(Hyp)[2]+x)}, aes(col="Hyperbolic"))
p+pHyp

PowStart <- nls2(PowFormula, data=xy, start=abc, algorithm="brute-force")
Pow <- nls2(PowFormula, data=xy, start=PowStart)
pPow <- stat_function(fun=function(x){coef(Pow)[1]+coef(Pow)[2]*x^coef(Pow)[3]}, aes(col="Power"))
p+pPow

source("3 nlsTab.R")

##########################################################################################################################################

response <- "funcMaxed70"
responseName <- "Number of functions above 70% threshold"

xy <- mfData[,c(paste(predictor), paste(response))]
colnames(xy) <- c("x", "y")
p <- qplot(x, y, data=xy) + labs(x=paste(eval(predictorName)), y=paste(eval(responseName)))
p

HoStart <- nls2(HoFormula, data=xy, start=a, algorithm="brute-force")
Ho <- nls2(HoFormula, data=xy, start=HoStart)
pHo <- stat_function(fun=function(x){coef(Ho)[1]+0*x}, aes(col="Null"))
p+pHo

LmStart <- nls2(LmFormula, data=xy, start=ab, algorithm="brute-force")
Lm <- nls2(LmFormula, data=xy, start=LmStart)
pLm <-stat_function(fun=function(x){coef(Lm)[1]+coef(Lm)[2]*x}, aes(col="Linear"))
p+pLm

LogStart <- nls2(LogFormula, data=xy, start=ab, algorithm="brute-force")
Log <- nls2(LogFormula, data=xy, start=LogStart)
pLog <- stat_function(fun=function(x){coef(Log)[1]+coef(Log)[2]*log(x+1)}, aes(col="Log"))
p+pLog

HypStart <- nls2(HypFormula, data=xy, start=ab, algorithm="brute-force")
Hyp <- nls2(HypFormula, data=xy, start=HypStart)
pHyp <- stat_function(fun=function(x){coef(Hyp)[1]*x/(coef(Hyp)[2]+x)}, aes(col="Hyperbolic"))
p+pHyp

PowStart <- nls2(PowFormula, data=xy, start=abc, algorithm="brute-force")
Pow <- nls2(PowFormula, data=xy, start=PowStart)
pPow <- stat_function(fun=function(x){coef(Pow)[1]+coef(Pow)[2]*x^coef(Pow)[3]}, aes(col="Power"))
p+pPow

source("3 nlsTab.R")

##########################################################################################################################################

response <- "funcMaxed80"
responseName <- "Number of functions above 80% threshold"

xy <- mfData[,c(paste(predictor), paste(response))]
colnames(xy) <- c("x", "y")
p <- qplot(x, y, data=xy) + labs(x=paste(eval(predictorName)), y=paste(eval(responseName)))
p

HoStart <- nls2(HoFormula, data=xy, start=a, algorithm="brute-force")
Ho <- nls2(HoFormula, data=xy, start=HoStart)
pHo <- stat_function(fun=function(x){coef(Ho)[1]+0*x}, aes(col="Null"))
p+pHo

LmStart <- nls2(LmFormula, data=xy, start=ab, algorithm="brute-force")
Lm <- nls2(LmFormula, data=xy, start=LmStart)
pLm <-stat_function(fun=function(x){coef(Lm)[1]+coef(Lm)[2]*x}, aes(col="Linear"))
p+pLm

LogStart <- nls2(LogFormula, data=xy, start=ab, algorithm="brute-force")
Log <- nls2(LogFormula, data=xy, start=LogStart)
pLog <- stat_function(fun=function(x){coef(Log)[1]+coef(Log)[2]*log(x+1)}, aes(col="Log"))
p+pLog

HypStart <- nls2(HypFormula, data=xy, start=ab, algorithm="brute-force")
Hyp <- nls2(HypFormula, data=xy, start=HypStart)
pHyp <- stat_function(fun=function(x){coef(Hyp)[1]*x/(coef(Hyp)[2]+x)}, aes(col="Hyperbolic"))
p+pHyp

PowStart <- nls2(PowFormula, data=xy, start=abc, algorithm="brute-force")
Pow <- nls2(PowFormula, data=xy, start=PowStart)
pPow <- stat_function(fun=function(x){coef(Pow)[1]+coef(Pow)[2]*x^coef(Pow)[3]}, aes(col="Power"))
p+pPow

source("3 nlsTab.R")

##########################################################################################################################################

response <- "funcMaxed90"
responseName <- "Number of functions above 90% threshold"

xy <- mfData[,c(paste(predictor), paste(response))]
colnames(xy) <- c("x", "y")
p <- qplot(x, y, data=xy) + labs(x=paste(eval(predictorName)), y=paste(eval(responseName)))
p

HoStart <- nls2(HoFormula, data=xy, start=a, algorithm="brute-force")
Ho <- nls2(HoFormula, data=xy, start=HoStart)
pHo <- stat_function(fun=function(x){coef(Ho)[1]+0*x}, aes(col="Null"))
p+pHo

LmStart <- nls2(LmFormula, data=xy, start=ab, algorithm="brute-force")
Lm <- nls2(LmFormula, data=xy, start=LmStart)
pLm <-stat_function(fun=function(x){coef(Lm)[1]+coef(Lm)[2]*x}, aes(col="Linear"))
p+pLm

LogStart <- nls2(LogFormula, data=xy, start=ab, algorithm="brute-force")
Log <- nls2(LogFormula, data=xy, start=LogStart)
pLog <- stat_function(fun=function(x){coef(Log)[1]+coef(Log)[2]*log(x+1)}, aes(col="Log"))
p+pLog

HypStart <- nls2(HypFormula, data=xy, start=ab, algorithm="brute-force")
Hyp <- nls2(HypFormula, data=xy, start=HypStart)
pHyp <- stat_function(fun=function(x){coef(Hyp)[1]*x/(coef(Hyp)[2]+x)}, aes(col="Hyperbolic"))
p+pHyp

PowStart <- nls2(PowFormula, data=xy, start=abc, algorithm="brute-force")
Pow <- nls2(PowFormula, data=xy, start=PowStart)
pPow <- stat_function(fun=function(x){coef(Pow)[1]+coef(Pow)[2]*x^coef(Pow)[3]}, aes(col="Power"))
p+pPow

source("3 nlsTab.R")

##########################################################################################################################################

# Write model selection table
write.csv(cum.nlsTab, "4 model selection table.csv")

##########################################################################################################################################

# Melt data for plotting
melt <- melt(mfData, id.vars=c(idVars, paste(predictor)))
melt$TrtPeg <- factor(melt$TrtPeg)

# Give pretty names for plotting
melt$variable2 <- melt$variable
melt$variable2 <- ifelse(melt$variable=="Epi", paste("paste('Epifauna abundance (# m'^-2, ')')"), paste(melt$variable2))
melt$variable2 <- ifelse(melt$variable=="EpiRich", paste("paste('Epifauna richness (taxa m'^-2, ')')"), paste(melt$variable2))
melt$variable2 <- ifelse(melt$variable=="Dsln", paste("paste('Dissolution (g d'^-1, ')')"), paste(melt$variable2))
melt$variable2 <- ifelse(melt$variable=="Sed", paste("paste('Sediment stabilization (', Delta, 'cm mo'^-1, ')')"), paste(melt$variable2))
melt$variable2 <- ifelse(melt$variable=="Nrsy", paste("paste('Nursery abundance (# m'^-2, ')')"), paste(melt$variable2))
melt$variable2 <- ifelse(melt$variable=="NrsyRich", paste("paste('Nursery richness (taxa m'^-2, ')')"), paste(melt$variable2))
melt$variable2 <- ifelse(melt$variable=="Dcmp", paste("paste('Decomposition (g mo'^-1, ')')"), paste(melt$variable2))
melt$variable2 <- ifelse(melt$variable=="meanFunction", paste("paste('Multifunctionality (%)')"), paste(melt$variable2))
melt$variable2 <- ifelse(melt$variable=="funcMaxed10", paste("paste('Number of functions'>= 'threshold')"), paste(melt$variable2))
melt$variable2 <- ifelse(melt$variable=="funcMaxed20", paste("paste('Number of functions'>= 'threshold')"), paste(melt$variable2))
melt$variable2 <- ifelse(melt$variable=="funcMaxed30", paste("paste('Number of functions'>= 'threshold')"), paste(melt$variable2))
melt$variable2 <- ifelse(melt$variable=="funcMaxed40", paste("paste('Number of functions'>= 'threshold')"), paste(melt$variable2))
melt$variable2 <- ifelse(melt$variable=="funcMaxed50", paste("paste('Number of functions'>= 'threshold')"), paste(melt$variable2))
melt$variable2 <- ifelse(melt$variable=="funcMaxed60", paste("paste('Number of functions'>= 'threshold')"), paste(melt$variable2))
melt$variable2 <- ifelse(melt$variable=="funcMaxed70", paste("paste('Number of functions'>= 'threshold')"), paste(melt$variable2))
melt$variable2 <- ifelse(melt$variable=="funcMaxed80", paste("paste('Number of functions'>= 'threshold')"), paste(melt$variable2))
melt$variable2 <- ifelse(melt$variable=="funcMaxed90", paste("paste('Number of functions'>= 'threshold')"), paste(melt$variable2))
melt$variable2 <- factor(melt$variable2, levels=c("paste('Epifauna abundance (# m'^-2, ')')", "paste('Epifauna richness (taxa m'^-2, ')')", "paste('Dissolution (g d'^-1, ')')", "paste('Sediment stabilization (', Delta, 'cm mo'^-1, ')')", "paste('Nursery abundance (# m'^-2, ')')", "paste('Nursery richness (taxa m'^-2, ')')", "paste('Decomposition (g mo'^-1, ')')", "paste('Multifunctionality (%)')", "paste('Number of functions'>= 'threshold')")) 

# Subset functions for which you plan to add lines
Epi <- subset(melt, melt$variable2=="paste('Epifauna abundance (# m'^-2, ')')")
Epi$variable2 <- factor(Epi$variable2)
EpiRich <- subset(melt, melt$variable2=="paste('Epifauna richness (taxa m'^-2, ')')")
EpiRich$variable2 <- factor(EpiRich$variable2)
Dsln <- subset(melt, melt$variable2=="paste('Dissolution (g d'^-1, ')')")
Dsln$variable2 <- factor(Dsln$variable2)
Sed <- subset(melt, melt$variable2=="paste('Sediment stabilization (', Delta, 'cm mo'^-1, ')')")
Sed$variable2 <- factor(Sed$variable2)
Nrsy <- subset(melt, melt$variable2=="paste('Nursery abundance (# m'^-2, ')')")
Nrsy$variable2 <- factor(Nrsy$variable2)
NrsyRich <- subset(melt, melt$variable2=="paste('Nursery richness (taxa m'^-2, ')')")
NrsyRich$variable2 <- factor(NrsyRich$variable2)
Dcmp <- subset(melt, melt$variable2=="paste('Decomposition (g mo'^-1, ')')")
Dcmp$variable2 <- factor(Dcmp$variable2)
meanFunction <- subset(melt, melt$variable2=="paste('Multifunctionality (%)')")
meanFunction$variable2 <- factor(meanFunction$variable2)

# Subset to plot points vs. threshold lines
funcs <- subset(melt, melt$variable2!="paste('Number of functions'>= 'threshold')")
funcs$variable2 <- factor(funcs$variable2)
tholds <- subset(melt, melt$variable2=="paste('Number of functions'>= 'threshold')")
tholds$variable2 <- factor(tholds$variable2)

# Subset to plot individual thresholds
tholds$thresholds <- as.character(tholds$variable)
tholds$thresholds <- as.numeric(with(tholds, gsub("funcMaxed", "", tholds$thresholds)))
t10 <- subset(tholds, thresholds == 10)
t20 <- subset(tholds, thresholds == 20)
t30 <- subset(tholds, thresholds == 30)
t40 <- subset(tholds, thresholds == 40)
t50 <- subset(tholds, thresholds == 50)
t60 <- subset(tholds, thresholds == 60)
t70 <- subset(tholds, thresholds == 70)
t80 <- subset(tholds, thresholds == 80)
t90 <- subset(tholds, thresholds == 90)

##########################################################################################################################################

# Generate plot
fig2 <- ggplot(aes(x=Gcvr, y=value), data=melt) + 
   facet_wrap(~variable2, scales="free", ncol=3, labeller=label_parsed) + 
   stat_function(data=Sed, fun=function(x){0*x}, size=0.25, lty=2, alpha=1) +
   theme_bw(base_size=6, base_family="Helvetica") +
   geom_point(data=funcs, aes(shape=factor(TrtPeg), fill=factor(TrtPeg)), stroke=0.25, color="black") + 
      scale_shape_manual(values=c(21, 25, 22, 21, 24, 23)) +  
      scale_fill_manual(values=c("red", "orange", "yellow", "green", "blue", "purple")) + 
      guides(shape=F, fill=F) +
   stat_function(data=Epi, fun=function(x){256.902*x/(x+1.037)}, size=0.25) + 
   stat_function(data=EpiRich, fun=function(x){6.535*x/(x+1.283)}, size=0.25) + 
   stat_function(data=Dsln, fun=function(x){9.352-0.313*log(x+1)}, size=0.25) + #DslnFlip?
   stat_function(data=Nrsy, fun=function(x){7.778*x/(x+37.253)}, size=0.25) + 
   stat_function(data=NrsyRich, fun=function(x){0.979+0.185*log(x+1)}, size=0.25) +
   stat_function(data=meanFunction, fun=function(x){0.364+0.082*log(x+1)}, size=0.25) + 
   stat_function(fun=function(x){4.144+0.690*log(x+1)}, data=t10, aes(colour=thresholds)) +
   stat_function(fun=function(x){3.769+0.728*log(x+1)}, data=t20, aes(colour=thresholds)) +
   stat_function(fun=function(x){3.435+0.751*log(x+1)}, data=t30, aes(colour=thresholds)) +
   stat_function(fun=function(x){3.498+0.683*log(x+1)}, data=t40, aes(colour=thresholds)) +
   stat_function(fun=function(x){2.380+1.672*x^0.188}, data=t50, aes(colour=thresholds)) +
   stat_function(fun=function(x){2.128+1.689*x^0.185}, data=t60, aes(colour=thresholds)) +
   stat_function(fun=function(x){1.907+0.769*log(x+1)}, data=t70, aes(colour=thresholds)) + 
   stat_function(fun=function(x){1.179+0.723*log(x+1)}, data=t80, aes(colour=thresholds)) +
   stat_function(fun=function(x){2.431*x/(x+4.868)}, data=t90, aes(colour=thresholds)) +
   xlab("\nGracilaria cover (%)") +
   theme(aspect.ratio=1,
      panel.background=element_rect(fill=NA),
      line=element_line(size=0.25),
      axis.title.y=element_blank(),
      axis.title.x=element_text(size=unit(6.5, "pt")),
      axis.ticks.length=unit(-0.05, "cm"),
      axis.text.x= element_text(margin = margin(t = 3)),
      axis.text.y=element_text(margin = margin (r=3)),
      panel.grid=element_blank(),
      strip.background=element_rect(fill=NA, colour=NA),
      axis.text=element_text(colour="black", size=unit(6, "pt")),
      strip.text=element_text(size=unit(6, "pt")),
      panel.spacing=unit(0.5, "lines"),
      panel.border=element_rect(fill=NA, color=NA), 
      axis.line=element_line(size=0.25),
      legend.justification="top",
      legend.position=c(0.9225, 0.08),
      legend.direction="horizontal",
      legend.background=element_rect(fill=NA),
      legend.margin=unit(0, "cm"),
      legend.text=element_text(size=unit(6, "pt")),
      legend.key=element_rect(fill = NA, colour = NA),
      legend.key.size = unit(1, "cm"),
      legend.key.height =  unit(0.15, "cm"),
      legend.key.width =   unit(0.2625, "cm"),
      legend.title=element_text(face="plain", size=unit(6, "pt"))) +
   scale_color_gradientn(name="Threshold (%)", 
      colours=rev(rainbow(5)), # 5
      guide=guide_colourbar(title.position="top"), breaks=c(10, 50, 90))

ggsave(fig2, file="5 figure.png", width=4.5, height=5.1, dpi=1200)
