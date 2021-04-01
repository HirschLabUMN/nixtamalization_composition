## Script written by Jonathan Renk
## 27 Mar 2020

## Clearing the global environment
rm(list=ls(all=TRUE))

## Loading in packages
library("ggplot2")
library("car")
library("lme4")
library("lmerTest")
library("dfoptim")

## Setting up the working directory
getwd()
setwd("/Users/jonathanrenk/Desktop/Discovery Work/Discovery Work")

## Loading in the data
data <- read.csv("data/master_compile_data_2.csv", header=T, stringsAsFactors=F)

## Setting variables as factors for genotype, subsample, and Rep
data[,1] <- as.factor(data[,1])
data[,2] <- as.factor(data[,2])
data[,3] <- as.factor(data[,3])

str(data)

### ALL OF THE DATA #### 
data <- read.csv("data/master_compile_data_2.csv", header=T, stringsAsFactors=F)

data[,1] <- as.factor(data[,1])
data[,2] <- as.factor(data[,2])
data[,3] <- as.factor(data[,3])

str(data)

traits <- colnames(data)[4:ncol(data)]

for(trait in traits){
  pdf(paste0("plots_all_", trait, ".pdf"), width = 15, height = 5)
  par(mfrow=c(1,2))
  
  # Histogram of raw data
  hist(data[,trait],
       main=paste0("Histogram of raw\n", trait, " values"),
       xlab=trait,
       col="cadetblue")
  
  # Stripchart of raw data across subsamples
  stripchart(data[,trait] ~data$Subsample,
             main=paste0("Stripchart of raw\n", trait, " values"),
             xlab="Subsample",
             ylab=trait,
             vertical=TRUE, 
             method="jitter",
             bg="cadetblue",
             pch=21
  )
  
  dev.off()
  
  if(paste0("stats_all_", trait, ".txt") %in% list.files()){
    system(paste0("rm stats_all_", trait, ".txt"))}
  # Summary statistics of the trait
  summary <- summary(data[,trait], )
  out <- capture.output(summary)
  cat(out, file=paste0("stats_all_",trait,".txt"), sep="\n", append=TRUE)
  
  # Test for normality
  normality <- shapiro.test(data[,trait])
  out <- capture.output(normality)
  cat(out, file=paste0("stats_all_",trait,".txt"), sep="\n", append=TRUE)
  
  # Run an ANOVA
  model <- lm(get(trait) ~ Genotype + Subsample + Rep + Genotype:Subsample + Genotype:Rep + Subsample:Rep, data=data)
  anova <- anova(model)
  out <- capture.output(anova)
  cat(out, file=paste0("stats_all_",trait,".txt"), sep="\n", append=TRUE)
  
  # Run a mixed model
  model.1 <- lmer(get(trait) ~ Genotype + Subsample + (1|Rep) + Genotype:Subsample+ (1|Genotype:Rep) + (1|Subsample:Rep), data = data, REML = TRUE)
  # Decreasing stopping tolerances
  strict_tol <- lmerControl(optCtrl=list(xtol_abs=1e-8, ftol_abs=1e-8))
  if (all(model.1@optinfo$optimizer=="nloptwrap")) {
    model <- update(model.1, control=strict_tol)
  }
  summary(model, correlation=FALSE)
  # ANOVA of fixed effects
  anova <- anova(model)
  out <- capture.output(anova)
  cat(out, file=paste0("stats_all_",trait,".txt"), sep="\n", append=TRUE)
  # Testing for random effects
  rand_eff <- rand(model)
  out <- capture.output(rand_eff)
  cat(out, file=paste0("stats_all_",trait,".txt"), sep="\n", append=TRUE)
  
  # Write out residuals from ANOVA
  write.table(resid(model), paste0("resids_all_", trait, ".csv"), col.names=F, row.names=F, sep=",")
  
  pdf(paste0("assumptions_all_", trait, ".pdf"), width = 15, height = 5)
  par(mfrow=c(1,3))
  
  # Model Fit with REML
  plot(fitted(model), residuals(model), pch=19, col="dark blue", ylab="Residuals", xlab="Predicted")
  abline(h=0,col="red", lwd=1, lty=1)
  # histogram of residuals
  hist(residuals(model),main="Histogram of residuals",freq=F, xlab="Residuals", ylab= "Freq", col="palegreen", col.main="darkblue")
  x=seq(-5e-15,9e-15,5e-15)
  curve(dnorm(x,mean(residuals(model)),sd(residuals(model))),add=T,lwd=2, col="red", lty=1)
  # qq plot
  qqPlot(residuals(model), pch=19, col="dark blue", col.lines="red", xlab="Pred quantiles", ylab="Obs quantiles") 
  
  dev.off()
  
}

### NIX ONLY ###

data <- read.csv("data/master_compile_data_2.csv", header=T, stringsAsFactors=F)

# subsetting the data
nix <- data[which(data$Subsample == "NIX"),]

nix[,1] <- as.factor(nix[,1])
nix[,2] <- as.factor(nix[,2])
nix[,3] <- as.factor(nix[,3])

str(nix)

traits <- colnames(nix)[4:ncol(nix)]

for(trait in traits){
  pdf(paste0("plots_nix", trait, ".pdf"), width = 15, height = 5)
  par(mfrow=c(1,2))
  
  # Histogram of raw data
  hist(nix[,trait],
       main=paste0("Histogram of raw\n", trait, " values"),
       xlab=trait,
       col="cadetblue")
  
  # Stripchart of raw data across subsamples
  stripchart(nix[,trait] ~nix$Subsample,
             main=paste0("Stripchart of raw\n", trait, " values"),
             xlab="Subsample",
             ylab=trait,
             vertical=TRUE, 
             method="jitter",
             bg="cadetblue",
             pch=21
  )
  
  dev.off()
  
  if(paste0("stats_nix", trait, ".txt") %in% list.files()){
    system(paste0("rm stats_nix", trait, ".txt"))}
  # Summary statistics of the trait
  summary <- summary(nix[,trait], )
  out <- capture.output(summary)
  cat(out, file=paste0("stats_nix",trait,".txt"), sep="\n", append=TRUE)
  
  # Test for normality
  normality <- shapiro.test(nix[,trait])
  out <- capture.output(normality)
  cat(out, file=paste0("stats_nix",trait,".txt"), sep="\n", append=TRUE)
  
  # Run an ANOVA
  model <- lm(get(trait) ~ Genotype + Rep, data=nix)
  anova <- anova(model)
  out <- capture.output(anova)
  cat(out, file=paste0("stats_nix",trait,".txt"), sep="\n", append=TRUE)
  
}

### WK ONLY ###
wk <- data[which(data$Subsample == "WK"),]

wk[,1] <- as.factor(wk[,1])
wk[,2] <- as.factor(wk[,2])
wk[,3] <- as.factor(wk[,3])

str(wk)

traits <- colnames(wk)[4:ncol(wk)]

for(trait in traits){
  pdf(paste0("plots_wk", trait, ".pdf"), width = 15, height = 5)
  par(mfrow=c(1,2))
  
  # Histogram of raw data
  hist(wk[,trait],
       main=paste0("Histogram of raw\n", trait, " values"),
       xlab=trait,
       col="cadetblue")
  
  # Stripchart of raw data across subsamples
  stripchart(wk[,trait] ~wk$Subsample,
             main=paste0("Stripchart of raw\n", trait, " values"),
             xlab="Subsample",
             ylab=trait,
             vertical=TRUE, 
             method="jitter",
             bg="cadetblue",
             pch=21
  )
  
  dev.off()
  
  if(paste0("stats_wk", trait, ".txt") %in% list.files()){
    system(paste0("rm stats_wk", trait, ".txt"))}
  # Summary statistics of the trait
  summary <- summary(wk[,trait], )
  out <- capture.output(summary)
  cat(out, file=paste0("stats_wk",trait,".txt"), sep="\n", append=TRUE)
  
  # Test for normality
  normality <- shapiro.test(wk[,trait])
  out <- capture.output(normality)
  cat(out, file=paste0("stats_wk",trait,".txt"), sep="\n", append=TRUE)
  
  # Run an ANOVA
  model <- lm(get(trait) ~ Genotype + Rep, data=wk)
  anova <- anova(model)
  out <- capture.output(anova)
  cat(out, file=paste0("stats_wk",trait,".txt"), sep="\n", append=TRUE)
  
}

### KRN ONLY ###
krn <- data[which(data$Subsample == "KRN"),]

krn[,1] <- as.factor(krn[,1])
krn[,2] <- as.factor(krn[,2])
krn[,3] <- as.factor(krn[,3])

str(krn)

traits <- colnames(krn)[4:ncol(krn)]

for(trait in traits){
  pdf(paste0("plots_krn", trait, ".pdf"), width = 15, height = 5)
  par(mfrow=c(1,2))
  
  # Histogram of raw data
  hist(krn[,trait],
       main=paste0("Histogram of raw\n", trait, " values"),
       xlab=trait,
       col="cadetblue")
  
  # Stripchart of raw data across subsamples
  stripchart(krn[,trait] ~krn$Subsample,
             main=paste0("Stripchart of raw\n", trait, " values"),
             xlab="Subsample",
             ylab=trait,
             vertical=TRUE, 
             method="jitter",
             bg="cadetblue",
             pch=21
  )
  
  dev.off()
  
  if(paste0("stats_krn", trait, ".txt") %in% list.files()){
    system(paste0("rm stats_krn", trait, ".txt"))}
  # Summary statistics of the trait
  summary <- summary(krn[,trait], )
  out <- capture.output(summary)
  cat(out, file=paste0("stats_krn",trait,".txt"), sep="\n", append=TRUE)
  
  # Test for normality
  normality <- shapiro.test(krn[,trait])
  out <- capture.output(normality)
  cat(out, file=paste0("stats_krn",trait,".txt"), sep="\n", append=TRUE)
  
  # Run an ANOVA
  model <- lm(get(trait) ~ Genotype + Rep, data=krn)
  anova <- anova(model)
  out <- capture.output(anova)
  cat(out, file=paste0("stats_krn",trait,".txt"), sep="\n", append=TRUE)
  
}

### Calculating emmeans ###
## Clearing the global environment
rm(list=ls(all=TRUE))

## Loading in packages
library("car")
library("lme4")
library("lmerTest")
library("dfoptim")
library("emmeans")
library("multcomp")

## Loading in the data
data <- read.csv("data/master_compile_data_2.csv", header=T, stringsAsFactors=F)

# removing nej and ws
data <- data[c(-721:-1200),]

# removing three genotypes for nitrogen and crude protein (B70, IBB14, PHK93)
#data <- data[c(-21,-38,-90,-141,-158,-210,-261,-278,-330,-381,-398,-450,-501,-518,-570,-621,-638,-690),]

data[,1] <- as.factor(data[,1])
data[,2] <- as.factor(data[,2])
data[,3] <- as.factor(data[,3])

str(data)
summary(data)

krn <- data[which(data$Subsample == "KRN"),]
krn[,1] <- as.factor(krn[,1])
krn[,2] <- as.factor(krn[,2])
krn[,3] <- as.factor(krn[,3])
str(krn)

nix <- data[which(data$Subsample == "NIX"),]
nix[,1] <- as.factor(nix[,1])
nix[,2] <- as.factor(nix[,2])
nix[,3] <- as.factor(nix[,3])
str(nix)

wk <- data[which(data$Subsample == "WK"),]
wk[,1] <- as.factor(wk[,1])
wk[,2] <- as.factor(wk[,2])
wk[,3] <- as.factor(wk[,3])
str(wk)

# switch for every trait and subsample of dataset #
model.1 <- lmer(Sucrose_g.100g_dwb ~ Genotype + (1|Rep), data = wk, REML = TRUE)
summary(model.1, correlation=FALSE)
anova(model.1) # fixed effects
rand(model.1)

# Decreasing stopping tolerances
strict_tol <- lmerControl(optCtrl=list(xtol_abs=1e-8, ftol_abs=1e-8))
if (all(model.1@optinfo$optimizer=="nloptwrap")) {
  model <- update(model.1, control=strict_tol)
}

# Emmeans
emm1 <- emmeans(model.1, ~ Genotype)
new_data <- as.data.frame(summary(emm1))[c('Genotype', 'emmean', 'SE', 'df', 'lower.CL', 'upper.CL')]
write.table(new_data, file="WK_Sucrose_dwb_emmeans.csv", col.names=F, row.names=F, sep=",")
marginal <- emmeans(model.1, ~ Subsample) # Running a Tukey's comparison
