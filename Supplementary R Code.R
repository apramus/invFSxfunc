# R code used to calculate multifunctionality, fit nonlinear models, and generated the figs presented in Ramus et al. 2017 PNAS
# Aaron P. Ramus (aaron.ramus@gmail.com)
# Last updated 20 Mar 2018

# Clear field
rm(list=ls(all.names=T))

setwd("~/Desktop/invFSxfunc-master")  # remove this at end

# Load required libraries
library(AICcmodavg)
library(ggplot2)
library(nls2)
library(broom)
library(multifunc)
library(ggpubr)
library(reshape2)
#library(cowplot)
#library(RColorBrewer)
#library(gridExtra)
#library(vegan)

# Read in data
meanPlots <- read.csv("mean plot-level responses.csv")
#mplots

# define id variables
id.vars <- c("Plot", "TrtPeg")

# select and name predictor variables
predictor <- "Gcvr"
predictor.name <- "Gracilaria cover (%)"

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

# create data frame for fits and plts
mfData <- meanPlots[,c(id.vars, eval(predictor), mfVars, "meanFunction")]

# Calculate multifunctionality thresholds (again, see Byrnes et al. 2014 Method Ecol Evol for usage of package 'multifunc')
mfThresh <- getFuncsMaxed(meanPlots, mfVars, threshmin=0.10, threshmax=0.90, threshstep=.1, maxN=9, prepend=id.vars)

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

# Source function modfit to extract information and parameters data frame to fit models
source("modfit.R")
source("parms.R")

# define models # formulas
HoFunc <- y ~ a
LmFunc <- y ~ a + b * x
LogFunc <- y ~ a + b * log(x + 1)
HypFunc <- y ~ a * x/(b + x)
PowFunc <- y ~ a + b * x^c

# create data frames to stuff things in later
cumStats <- NULL
cumFits <- NULL

# select and name individual response variable to be analyzed
response <- "Epi"
response.name <- "Epifauna abundance (# m^-2)"

# create data frame for analysis
xy <- mfData[,c(paste(predictor), paste(response))]
colnames(xy) <- c("x", "y")
p <- qplot(x, y, data=xy) + labs(x=paste(eval(predictor.name)), y=paste(eval(response.name)))
p

# fit models
HoStart <- nls2(HoFunc, data=xy, start=a, algorithm="brute-force")
Ho <- nls2(HoFunc, data=xy, start=HoStart)
pHo <- stat_function(fun=function(x){coef(Ho)[1]+0*x}, aes(col="Null"))
p+pHo

LmStart <- nls2(LmFunc, data=xy, start=ab, algorithm="brute-force")
Lm <- nls2(LmFunc, data=xy, start=LmStart)
pLm <-stat_function(fun=function(x){coef(Lm)[1]+coef(Lm)[2]*x}, aes(col="Linear"))
p+pLm

LogStart <- nls2(LogFunc, data=xy, start=ab, algorithm="brute-force")
Log <- nls2(LogFunc, data=xy, start=LogStart)
pLog <- stat_function(fun=function(x){coef(Log)[1]+coef(Log)[2]*log(x+1)}, aes(col="Log"))
p+pLog

HypStart <- nls2(HypFunc, data=xy, start=ab, algorithm="brute-force")
Hyp <- nls2(HypFunc, data=xy, start=HypStart)
#Hyp <- nls2(HypFunc, data=xy, start=c(a=75, b=1)) #, nls.control(maxiter=200))
pHyp <- stat_function(fun=function(x){coef(Hyp)[1]*x/(coef(Hyp)[2]+x)}, aes(col="Hyperbolic"))
p+pHyp

PowStart <- nls2(PowFunc, data=xy, start=abc, algorithm="brute-force")
Pow <- nls2(PowFunc, data=xy, start=PowStart)
#Pow <- nls2(PowFunc, data=xy, start=c(a=0.2, b=1, c=1))
pPow <- stat_function(fun=function(x){coef(Pow)[1]+coef(Pow)[2]*x^coef(Pow)[3]}, aes(col="Power"))
p+pPow

source("modComps.R")

write.csv(cumStats, "Model Selection Table.csv")





# Melt for plotting

melt <- melt(mfData, id.vars=c(id.vars, paste(predictor)))

str(melt)
melt$TrtPeg <- factor(melt$TrtPeg)

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

melt$variable2 <- factor(melt$variable2)
levels(melt$variable2)

melt$variable2 <- factor(melt$variable2, levels=c("paste('Epifauna abundance (# m'^-2, ')')", "paste('Epifauna richness (taxa m'^-2, ')')", "paste('Dissolution (g d'^-1, ')')", "paste('Sediment stabilization (', Delta, 'cm mo'^-1, ')')", "paste('Nursery abundance (# m'^-2, ')')", "paste('Nursery richness (taxa m'^-2, ')')", "paste('Decomposition (g mo'^-1, ')')", "paste('Multifunctionality (%)')", "paste('Number of functions'>= 'threshold')")) 

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


funcs <- subset(melt, melt$variable2!="paste('Number of functions'>= 'threshold')")
funcs$variable2 <- factor(funcs$variable2)

tholds <- subset(melt, melt$variable2=="paste('Number of functions'>= 'threshold')")
tholds$variable2 <- factor(tholds$variable2)

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

shapes=c(21, 25, 22, 21, 24, 23)
cols=c("red", "orange", "yellow", "green", "blue", "purple")

# Create Fig2
Fig2 <- ggplot(aes(x=Gcvr, y=value), data=melt) + 
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
   stat_function(fun=function(x){4.14+0.69*log(x+1)}, data=t10, aes(colour=thresholds)) +
   stat_function(fun=function(x){3.77+0.73*log(x+1)}, data=t20, aes(colour=thresholds)) +
   stat_function(fun=function(x){3.44+0.75*log(x+1)}, data=t30, aes(colour=thresholds)) +
   stat_function(fun=function(x){3.50+0.68*log(x+1)}, data=t40, aes(colour=thresholds)) +
   stat_function(fun=function(x){2.38+1.67*x^0.19}, data=t50, aes(colour=thresholds)) +
   stat_function(fun=function(x){2.13+1.69*x^0.19}, data=t60, aes(colour=thresholds)) +
   stat_function(fun=function(x){1.91+0.77*log(x+1)}, data=t70, aes(colour=thresholds)) + 
   stat_function(fun=function(x){1.18+0.72*log(x+1)}, data=t80, aes(colour=thresholds)) +
   stat_function(fun=function(x){2.43*x/(x+4.87)}, data=t90, aes(colour=thresholds)) +
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
      #legend.margin=unit(0, "cm"),
      legend.text=element_text(size=unit(6, "pt")),
      legend.key=element_rect(fill = NA, colour = NA),
      #legend.key.size = unit(1, "cm"),
      legend.key.height =  unit(0.15, "cm"),
      legend.key.width =   unit(0.2625, "cm"),
      legend.title=element_text(face="plain", size=unit(6, "pt"))) +
      scale_color_gradientn(name="Threshold (%)", 
         colours=rev(rainbow(5)), # 5
         guide=guide_colourbar(title.position="top"), breaks=c(10, 50, 90))
Fig2

ggsave(Fig2, file="Fig. 2. Density-dependent impacts of an invasive foundation species on coastal ecosystem functions.pdf", width=4.5, height=5.1, dpi=2400)
   

# Create Fig2
Fig2 <- ggplot(aes(x=Gcvr, y=value), data=melt) + 
   facet_wrap(~variable2, scales="free", ncol=3, labeller=label_parsed) + 
   stat_function(data=Sed, fun=function(x){0*x}, size=0.25, lty=2, alpha=1) +
   theme_pubr(base_size=6) + 
   geom_point(data=funcs, aes(shape=factor(TrtPeg), fill=factor(TrtPeg)), stroke=0.25, color="black") + 
      scale_shape_manual(values=c(16, 25, 22, 21, 24, 23)) +  
      scale_fill_manual(values=c("red", "orange", "yellow", "green", "blue", "purple")) + 
      guides(shape=F, fill=F) +
   stat_function(data=Epi, fun=function(x){256.902*x/(x+1.037)}, size=0.25) + 
   stat_function(data=EpiRich, fun=function(x){6.535*x/(x+1.283)}, size=0.25) + 
   stat_function(data=Dsln, fun=function(x){9.352-0.313*log(x+1)}, size=0.25) + #DslnFlip?
   stat_function(data=Nrsy, fun=function(x){7.778*x/(x+37.253)}, size=0.25) + 
   stat_function(data=NrsyRich, fun=function(x){0.979+0.185*log(x+1)}, size=0.25) +
   stat_function(data=meanFunction, fun=function(x){0.364+0.082*log(x+1)}, size=0.25) + 
   stat_function(fun=function(x){4.14+0.69*log(x+1)}, data=t10, aes(colour=thresholds)) +
   stat_function(fun=function(x){3.77+0.73*log(x+1)}, data=t20, aes(colour=thresholds)) +
   stat_function(fun=function(x){3.44+0.75*log(x+1)}, data=t30, aes(colour=thresholds)) +
   stat_function(fun=function(x){3.50+0.68*log(x+1)}, data=t40, aes(colour=thresholds)) +
   stat_function(fun=function(x){2.38+1.67*x^0.19}, data=t50, aes(colour=thresholds)) +
   stat_function(fun=function(x){2.13+1.69*x^0.19}, data=t60, aes(colour=thresholds)) +
   stat_function(fun=function(x){1.91+0.77*log(x+1)}, data=t70, aes(colour=thresholds)) + 
   stat_function(fun=function(x){1.18+0.72*log(x+1)}, data=t80, aes(colour=thresholds)) +
   stat_function(fun=function(x){2.43*x/(x+4.87)}, data=t90, aes(colour=thresholds)) +
   xlab("\nGracilaria cover (%)") +
      theme(aspect.ratio=1, 
         strip.background=element_rect(fill=NA, colour=NA),
         legend.justification="top",
         legend.position=c(0.9225, 0.08),
         legend.direction="horizontal",
         legend.background=element_rect(fill=NA),
         #legend.margin=unit(0, "cm"),
         legend.text=element_text(size=unit(6, "pt")),
         legend.key=element_rect(fill = NA, colour = NA),
         #legend.key.size = unit(1, "cm"),
         legend.key.height =  unit(0.15, "cm"),
         legend.key.width =   unit(0.2625, "cm"),
         legend.title=element_text(face="plain", size=unit(6, "pt"))) +
      scale_color_gradientn(
         name="Threshold (%)", 
         colours=rev(rainbow(5)), # 5
         guide=guide_colourbar(title.position="top"), 
         breaks=c(10, 50, 90)
         )
Fig2

ggsave(Fig2, file="Fig. 2. Density-dependent impacts of an invasive foundation species on coastal ecosystem functions.pdf", width=4.5, height=5.1, dpi=2400)
   
