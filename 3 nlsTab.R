##########################################################################################################################################
#                                                                                                                                        #
# Code to generate a model selection table for fitted nls models.                                                                        #
#                                                                                                                                        #
# This code works as a 'manual loop' for lack of a better description. Much of it is dedicated to combining and reorganizing             #
# infomration, although there are a few calculations.                                                                                    #
#                                                                                                                                        #
# Developed by Aaron Ramus (aaron.ramus@gmail.com). Last updated 21 March 2018.                                                          #
#                                                                                                                                        #
##########################################################################################################################################

# Create 'nlsTab' from aictab (package AICcmodavg)
nlsTab <- as.data.frame(aictab(list(Ho, Lm, Log, Hyp, Pow), c("Null", "Linear", "Log", "Hyperbolic", "Power")))

# Reorganize nlsTab and formula names
colnames(nlsTab)[1] <- "Model"
nlsTab$Formula <- as.character(nlsTab$Model)
nlsTab$Formula <- ifelse(nlsTab$Formula=="Null", paste("y ~ a"), nlsTab$Formula)
nlsTab$Formula <- ifelse(nlsTab$Formula=="Linear", paste("y ~ a + b * x"), nlsTab$Formula)
nlsTab$Formula <- ifelse(nlsTab$Formula=="Log", paste("y ~ a + b * log(x + 1)"), nlsTab$Formula)
nlsTab$Formula <- ifelse(nlsTab$Formula=="Hyperbolic", paste("y ~ a * x/(b + x)"), nlsTab$Formula)
nlsTab$Formula <- ifelse(nlsTab$Formula=="Power", paste("y ~ a + b * x^c"), nlsTab$Formula)
nlsTab$Formula <- factor(nlsTab$Formula)

# Return information about fitted nls models
HoFit <- getModFit(Ho)
LmFit<- getModFit(Lm)
LogFit <- getModFit(Log)
HypFit <- getModFit(Hyp)
PowFit <- getModFit(Pow)

# Bind returned information into a data frame, rename columns, and merge with nlsTab
modFits <- as.data.frame(rbind(HoFit, LmFit, LogFit, HypFit, PowFit))
colnames(modFits)[c(1,8,9,15)] <- c("Formula", "RSE", "df", "RSS")
nlsTab <- merge(nlsTab, modFits)

# Perform calculations, rename columns, and order table
nlsTab$nullModelLik <- nlsTab$ModelLik[1]
nlsTab$Response.variable <- eval(responseName)
nlsTab$Formula <- with(nlsTab, gsub("y ~ ", "", nlsTab$Formula))
nlsTab$Formula <- with(nlsTab, gsub("x", eval(predictor), nlsTab$Formula))
nlsTab$N <- nrow(xy)
nlsTab$R2 <- 1-(nlsTab$RSS/nlsTab$TSS)
nlsTab$adjR2 <- 1-((nlsTab$RSS/nlsTab$df)/(nlsTab$TSS/(nlsTab$N-1)))
nlsTab$genR2 <- 1-(nlsTab$nullModelLik/nlsTab$ModelLik)^(2/nlsTab$N)
nlsTab <- nlsTab[order(nlsTab$Delta_AICc),]

# Select columns to display in table, then round them and make pretty
nlsTab <- nlsTab[,c("Response.variable", "Model", "Formula", "K",  "a", "b", "c", "AICc", "Delta_AICc", "AICcWt", "LL", "df", "fvalue", "pvalue", "adjR2")]
colnames(nlsTab) <- c("Response.variable", "Model", "Formula", "K", "a", "b", "c", "AICc", "deltaAICc", "AICcWt", "LL", "df", "F", "P", "adjR2")
nlsTab$P <- ifelse(nlsTab$Model=="Null", 1, nlsTab$P)
nlsTab$'F' <- ifelse(nlsTab$Model=="Null", 0, nlsTab$'F')
nlsTab$PInitial <- nlsTab$P
nlsTab$P <- round(nlsTab$P, 3)
nlsTab$P <- ifelse(is.na(nlsTab$P), 1, nlsTab$P)
nlsTab$P <- ifelse(nlsTab$P==0, paste("<0.001"), nlsTab$P)
nlsTab[,c("a", "b", "c", "AICc", "deltaAICc", "AICcWt", "LL", "F",  "adjR2")] <- apply(nlsTab[,c("a", "b", "c", "AICc", "deltaAICc", "AICcWt", "LL", "F",  "adjR2")], 2, function(x) {round(x, 3)})
nlsTab <- nlsTab[,c("Response.variable", "Model", "Formula", "K", "a", "b", "c", "AICc", "deltaAICc", "AICcWt", "LL", "df", "F", "P", "adjR2")]

# Cumulatively bind to create full model selection table
cum.nlsTab <- as.data.frame(rbind(cum.nlsTab, nlsTab))

# Plot fitted models
pAll <- p + pHo + pLm + pLog + pHyp + pPow + scale_color_manual("Model", breaks=c("Null", "Linear", "Log", "Hyperbolic", "Power"), limits=c("Null", "Linear", "Log", "Hyperbolic", "Power"), values=rainbow(7)) + ylim(min(xy$y), max(xy$y))

# Display things
print(pAll)
print(nlsTab)

# Remove objects before fitting next function/response variable
rm(response, responseName, xy, p, HoStart, Ho, pHo, LmStart, Lm, pLm, LogStart, Log, pLog, HypStart, Hyp, pHyp, PowStart, Pow, pPow, nlsTab, HoFit, LmFit, LogFit, HypFit, PowFit, modFits, pAll)
