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

### ALL OF THE DATA ### only for HPLC traits
traits <- colnames(data)[9:16]

for(trait in traits){
  pdf(paste0("plots_hplc_", trait, ".pdf"), width = 15, height = 5)
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
  
  if(paste0("stats_hplc_", trait, ".txt") %in% list.files()){
    system(paste0("rm stats_hplc_", trait, ".txt"))}
  # Summary statistics of the trait
  summary <- summary(data[,trait], )
  out <- capture.output(summary)
  cat(out, file=paste0("stats_hplc_",trait,".txt"), sep="\n", append=TRUE)
  
  # Test for normality
  normality <- shapiro.test(data[,trait])
  out <- capture.output(normality)
  cat(out, file=paste0("stats_hplc_",trait,".txt"), sep="\n", append=TRUE)
  
  # Run an ANOVA
  model <- lm(get(trait) ~ Genotype + Subsample + Rep + Genotype:Subsample + Genotype:Rep + Subsample:Rep, data=data)
  anova <- anova(model)
  out <- capture.output(anova)
  cat(out, file=paste0("stats_hplc_",trait,".txt"), sep="\n", append=TRUE)
  
  # Run a mixed model
  model.1 <- lmer(get(trait) ~ Genotype + Subsample + (1|Rep) + Genotype:Subsample, data = data, REML = TRUE)
  # Decreasing stopping tolerances
  strict_tol <- lmerControl(optCtrl=list(xtol_abs=1e-8, ftol_abs=1e-8))
  if (all(model.1@optinfo$optimizer=="nloptwrap")) {
    model <- update(model.1, control=strict_tol)
  }
  summary(model, correlation=FALSE)
  # ANOVA of fixed effects
  anova <- anova(model)
  out <- capture.output(anova)
  cat(out, file=paste0("stats_hplc_",trait,".txt"), sep="\n", append=TRUE)
  # Testing for random effects
  rand_eff <- rand(model)
  out <- capture.output(rand_eff)
  cat(out, file=paste0("stats_hplc_",trait,".txt"), sep="\n", append=TRUE)
  
  # Write out residuals from ANOVA
  write.table(resid(model), paste0("resids_hplc_", trait, ".csv"), col.names=F, row.names=F, sep=",")
  
  pdf(paste0("assumptions_hplc_", trait, ".pdf"), width = 15, height = 5)
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

### ALL OF THE DATA excluding solutions ##### Loading in the data
data <- read.csv("data/master_compile_data_2.csv", header=T, stringsAsFactors=F)

# removing nej and ws
data <- data[c(-721:-1200),]

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

### NEJ ONLY ###
nej <- data[which(data$Subsample == "NEJ"),]

# Delete traits with no values
nej <- nej[,c(-4:-8)]

nej[,1] <- as.factor(nej[,1])
nej[,2] <- as.factor(nej[,2])
nej[,3] <- as.factor(nej[,3])

str(nej)

traits <- colnames(nej)[4:ncol(nej)]

for(trait in traits){
  pdf(paste0("plots_nej", trait, ".pdf"), width = 15, height = 5)
  par(mfrow=c(1,2))
  
  # Histogram of raw data
  hist(nej[,trait],
       main=paste0("Histogram of raw\n", trait, " values"),
       xlab=trait,
       col="cadetblue")
  
  # Stripchart of raw data across subsamples
  stripchart(nej[,trait] ~nej$Subsample,
             main=paste0("Stripchart of raw\n", trait, " values"),
             xlab="Subsample",
             ylab=trait,
             vertical=TRUE, 
             method="jitter",
             bg="cadetblue",
             pch=21
  )
  
  dev.off()
  
  if(paste0("stats_nej", trait, ".txt") %in% list.files()){
    system(paste0("rm stats_nej", trait, ".txt"))}
  # Summary statistics of the trait
  summary <- summary(nej[,trait], )
  out <- capture.output(summary)
  cat(out, file=paste0("stats_nej",trait,".txt"), sep="\n", append=TRUE)
  
  # Test for normality
  normality <- shapiro.test(nej[,trait])
  out <- capture.output(normality)
  cat(out, file=paste0("stats_nej",trait,".txt"), sep="\n", append=TRUE)
  
  # Run an ANOVA
  model <- lm(get(trait) ~ Genotype + Rep, data=nej)
  anova <- anova(model)
  out <- capture.output(anova)
  cat(out, file=paste0("stats_nej",trait,".txt"), sep="\n", append=TRUE)
  
}

### WS ONLY ###
ws <- data[which(data$Subsample == "WS"),]

# Delete traits with no values
ws <- ws[,c(-4:-8)]

ws[,1] <- as.factor(ws[,1])
ws[,2] <- as.factor(ws[,2])
ws[,3] <- as.factor(ws[,3])

str(ws)

traits <- colnames(ws)[4:ncol(ws)]

for(trait in traits){
  pdf(paste0("plots_ws", trait, ".pdf"), width = 15, height = 5)
  par(mfrow=c(1,2))
  
  # Histogram of raw data
  hist(ws[,trait],
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
  
  if(paste0("stats_ws", trait, ".txt") %in% list.files()){
    system(paste0("rm stats_ws", trait, ".txt"))}
  # Summary statistics of the trait
  summary <- summary(ws[,trait], )
  out <- capture.output(summary)
  cat(out, file=paste0("stats_ws",trait,".txt"), sep="\n", append=TRUE)
  
  # Test for normality
  normality <- shapiro.test(ws[,trait])
  out <- capture.output(normality)
  cat(out, file=paste0("stats_ws",trait,".txt"), sep="\n", append=TRUE)
  
  # Run an ANOVA
  model <- lm(get(trait) ~ Genotype + Rep, data=ws)
  anova <- anova(model)
  out <- capture.output(anova)
  cat(out, file=paste0("stats_ws",trait,".txt"), sep="\n", append=TRUE)
  
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
