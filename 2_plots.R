### Script written by Jonathan Renk
### 1 Apr 2021

### Figure 2a ###
## Clearing the global environment
rm(list=ls(all=TRUE))

## Loading in packages
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(agricolae)

### Extreme Selection of Inbreds ###
# loading in data 
data <- read.csv("data/selects_figure_2.csv", header = TRUE)
str(data)

# convert to long format
datal <- melt(data, id.vars="Extreme")
print(datal)

### Violin plots of extremes selection ###
plot <- ggplot(data = datal, aes(x=Extreme, y=value)) + 
  geom_violin() +
  geom_boxplot(aes(group=Extreme), width=0.1, color="black", position = position_dodge(width =0.9)) +
  labs(y = "Value (%)", x = "") +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 6)) +
  facet_wrap( ~ variable, scales="free", labeller = labeller(variable = c("Protein.As.is" = "Protein As Is","Starch.As.is" = "Starch As Is","Total_Sugars" = "Total Sugars"))) +
  theme(legend.position="none")

ggsave(plot, filename = "figure2_v2.pdf", width = 3.5, height = 3, units = c("in"))

### Figure 2b-d ###
## Clearing the global environment
rm(list=ls(all=TRUE))

## Setting up the working directory
data <- read.csv("data/nir_vs_lab.csv", header=T, stringsAsFactors=F)

#removing weird columns and rows
data <- data[c(-121),]
data <- data[,c(-8:-13)]

## Loading in packages
library(ggplot2)

## Loading in the data
#protein
str(data)

p <- ggplot(data = data, aes(x = Protein_Lab, y = Protein_NIR)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point()+
  ggtitle("Protein") +
  xlab("Lab (%)") +
  ylab("NIR (%)") +
  theme_classic()
p

lm_eqn <- function(data){
  m <- lm(Protein_NIR ~ Protein_Lab, data);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(sqrt(summary(m)$r.squared), digits = 3)))
  as.character(as.expression(eq));
}

p1 <- p + annotate("text", x = 9.5, y = 15, label = lm_eqn(data), parse = TRUE, size = 3)
p1
ggsave(p1, filename = "protein_nir.pdf", width = 3.5, height = 3, units = c("in"))

#starch
p <- ggplot(data = data, aes(x = Starch_Lab, y = Starch_NIR)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point()+
  ggtitle("Starch") +
  xlab("Lab (%)") +
  ylab("NIR (%)") +
  theme_classic()
p

lm_eqn <- function(data){
  m <- lm(Starch_NIR ~ Starch_Lab, data);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(sqrt(summary(m)$r.squared), digits = 3)))
  as.character(as.expression(eq));
}

p1 <- p + annotate("text", x = 64, y = 74, label = lm_eqn(data), parse = TRUE, size = 3)
p1
ggsave(p1, filename = "starch_nir.pdf", width = 3.5, height = 3.0, units = c("in"))

#total sugars
p <- ggplot(data = data, aes(x = Total_Sugars_Lab, y = Total_Sugars_NIR)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point()+
  ggtitle("Total Sugars") +
  xlab("Lab (g/100g flour)") +
  ylab("NIR (%)") +
  theme_classic()
p

lm_eqn <- function(data){
  m <- lm(Total_Sugars_NIR ~ Total_Sugars_Lab, data);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(sqrt(summary(m)$r.squared), digits = 3)))
  as.character(as.expression(eq));
}

p1 <- p + annotate("text", x = 27, y = 4.3, label = lm_eqn(data), parse = TRUE, size = 3)
p1
ggsave(p1, filename = "total_sugars_nir.pdf", width = 3.5, height = 3.0, units = c("in"))

### Figure 3 ###
## Clearing the global environment
rm(list=ls(all=TRUE))

## Loading in packages
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(agricolae)

## Boxplots of All of the data (excluding solutions)
data.1 <- read.csv("data/all_emmeans_final.csv", header = T, stringsAsFactors = F)
str(data.1)
#removing nej and ws
data.1 <- data.1[c(-361:-600),]
#setting subsample as factor 
data.1[,2] <- as.factor(data.1[,2])
data.1[,2] <- factor(data.1[,2], levels = c("RK", "CK", "SK"))
data.1 <- data.1[,c(-1,-3,-6,-8,-10:-12,-14,-16:-17)]

my_data <- data.1[, c(1, 3, 2, 6, 4, 5, 7)]

# convert to long format
datal.1 <- melt(my_data, id.vars="Subsample")
print(datal.1)

ggplot(data = datal.1, aes(x=variable, y=value, na.rm = TRUE)) + 
  geom_violin(aes(fill=Subsample), na.rm = TRUE, alpha = .75) +
  geom_boxplot(aes(group=Subsample), width=0.1, color="black", position = position_dodge(width =0.9), na.rm = TRUE) +
  facet_wrap( ~ variable, scales="free", labeller = labeller(variable = c("Crude.Protein" = "Crude Protein","Asparagine.dwb" = "Asparagine", "Starch.dwb" = "Starch", "Fructose.dwb" = "Fructose", "Glucose.dwb" = "Glucose", "Sucrose.dwb" = "Sucrose"))) +
  labs(y = "Value", x = "") +
  scale_fill_manual(values=c("gold", "darkred", "lightblue")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(legend.position="bottom")

testMeans <- function(datal.1){
  
  # Fit model to data for a single trait
  model <- lm(value ~ Subsample, data=datal.1)
  
  # Perform Tukey HSD
  test <- HSD.test(model, "Subsample")
  
  # Grab the significance letters from the tukey results
  datal.1$SigGroup <- test$groups$groups[datal.1$Subsample]
  out <- test$groups
  out$Subsample <- factor(rownames(out), levels=c("RK", "CK", "SK"))
  out <- out[order(out$Subsample), ]
  out$Trait <- unique(datal.1$variable)
  out$Mean <- out$value
  out$value <- max(datal.1$value)*1.15 
  
  # Calculate % change in phenotypes relative to your base (KRN in this case)
  BNIX <- (out["CK", "Mean"] - out["RK", "Mean"]) / out["RK", "Mean"]
  BWK <- (out["SK", "Mean"] - out["RK", "Mean"]) / out["RK", "Mean"]
  
  # Paste together the % change and the significance letters.  The empty "" is to hold space above KRN
  out$Label <- paste(c("", paste0(c(round(BNIX*100), round(BWK*100)), "%")),
                     out$groups,
                     sep="\n")
  out
}

# Use dplyr to apply the testMeans() function to each trait individually
diffs <- ddply(datal.1, .(variable), testMeans)

# Replacing value for Crude protein NAs 
diffs[1,2] = 14.46875*1.15 #KRN Crude Protein
diffs[2,2] = 14.46875*1.15 #NIX Crude Protein
diffs[3,2] = 14.46875*1.15 #WK Crude Protein

# Plot the figure
plot <- ggplot(datal.1, aes(Subsample, value, fill=Subsample)) +
  geom_violin(alpha=0.75) +
  stat_summary(fun.y="median", geom="point", pch='_', size=10, color="black", fill="white", na.rm = TRUE) +
  geom_boxplot(aes(group=Subsample), width=0.1, color="black", fill="white", position = position_dodge(width =0.9), na.rm = TRUE) +
  #geom_jitter(aes(color=Subsample), size=1.25, fill="white", pch=21, width=0.025, height=0) +
  facet_wrap(~variable, scales="free_y", ncol=6, labeller = labeller(variable = c("Crude.Protein" = "Crude Protein","Asparagine.dwb" = "Asparagine", "Starch.dwb" = "Starch", "Fructose.dwb" = "Fructose", "Glucose.dwb" = "Glucose", "Sucrose.dwb" = "Sucrose"))) +
  scale_fill_manual(values=c("gold", "darkred", "lightblue")) +
  scale_color_manual(values=c("gold", "darkred", "lightblue")) +
  geom_text(aes(label=Label), dat=diffs, color="black", vjust=1, size=2) + # ** This adds the labels**
  labs(y="Value (%)", x="") +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(legend.position="bottom")

ggsave(plot, filename = "figure3_v3.pdf", width = 7, height = 3, units = c("in"))

### Figure 4 ###
## Clearing the global environment
rm(list=ls(all=TRUE))

## Loading in packages
library(ggplot2)
library(reshape2)

par(mfrow=c(2,3), mai = c(0.25, 0.5, 0.25, 0.1))     

## Loading in the emmeans dataset
data.1 <- read.csv("data/all_emmeans_final.csv", header = T, stringsAsFactors = F)
#removing nej and ws
data.1 <- data.1[c(-361:-600),]
str(data.1)

## Setting variables as factors for genotype and subsample type
data.1[,1] <- as.factor(data.1[,1])
data.1[,2] <- as.factor(data.1[,2])
data.1[,2] <- factor(data.1[,2], levels = c("RK", "CK", "SK"))
str(data.1)

## Regression of all the data
model <- lm(Glucose.dwb ~ Sample, data = data.1)
summary(model)

model2 <- lm(Fructose.dwb ~ Sample, data = data.1)
summary(model2)

model3 <- lm(Sucrose.dwb ~ Sample, data = data.1)
summary(model3)

model4 <- lm(Asparagine.dwb ~ Sample, data = data.1)
summary(model4)

model5 <- lm(Crude.Protein ~ Sample, data = data.1)
summary(model5)

model6 <- lm(Starch.dwb ~ Sample, data = data.1)
summary(model6)

## Regression on all of the genotypes for glucose
genotypes <- unique(data.1$Genotype)
reg <- matrix(nrow = length(genotypes), ncol = 5, dimnames = list(genotypes, c("beta0", "beta1", "MSE", "Rsq", "pval")))

for(i in 1:length(genotypes)){
  temp <- data.1[data.1$Genotype == genotypes[i],]
  glucose_regression <- lm(temp[,9] ~ temp[,16])
  reg[i, "beta0"] <- glucose_regression$coefficients[1]
  reg[i, "beta1"] <- glucose_regression$coefficients[2]
  reg[i, "MSE"] <- sum(glucose_regression$residuals^2) / glucose_regression$df.residual
  reg[i, "Rsq"] <- summary(glucose_regression)$r.squared
  reg[i, "pval"] <- summary(glucose_regression)$coefficients[2,4]
}  

spread <- range(data.1[,9])
stripchart(data.1$Glucose.dwb ~ data.1$Subsample, 
           vertical = TRUE, 
           main = "",
           xlab="",
           ylab="",
           type="n")
for(i in 1:nrow(reg)){
  abline(a=reg[i,1], b=(reg[i,2]))
}
abline(a=model$coefficients[1], b=model$coefficients[2], col="red", lwd=3)
text(x=2.25, y=20,"p-value: <2e-16 ***")

## Regression on all samples for fructose
reg <- matrix(nrow = length(genotypes), ncol = 5, dimnames = list(genotypes, c("beta0", "beta1", "MSE", "Rsq", "pval")))

for(i in 1:length(genotypes)){
  temp <- data.1[data.1$Genotype == genotypes[i],]
  fructose_regression <- lm(temp[,7] ~ temp[,16])
  reg[i, "beta0"] <- fructose_regression$coefficients[1]
  reg[i, "beta1"] <- fructose_regression$coefficients[2]
  reg[i, "MSE"] <- sum(fructose_regression$residuals^2) / fructose_regression$df.residual
  reg[i, "Rsq"] <- summary(fructose_regression)$r.squared
  reg[i, "pval"] <- summary(fructose_regression)$coefficients[2,4]
}  

spread <- range(data.1[,7])
stripchart(data.1$Fructose.dwb ~ data.1$Subsample, 
           vertical = TRUE, 
           main = "",
           xlab="",
           ylab="",
           type="n")
for(i in 1:nrow(reg)){
  abline(a=reg[i,1], b=(reg[i,2]))
}
abline(a=model2$coefficients[1], b=model2$coefficients[2], col="red", lwd=3)
text(x=2.25, y=10,"p-value: <2e-16 ***")

## Regression on all samples for sucrose
reg <- matrix(nrow = length(genotypes), ncol = 5, dimnames = list(genotypes, c("beta0", "beta1", "MSE", "Rsq", "pval")))

for(i in 1:length(genotypes)){
  temp <- data.1[data.1$Genotype == genotypes[i],]
  sucrose_regression <- lm(temp[,15] ~ temp[,16])
  reg[i, "beta0"] <- sucrose_regression$coefficients[1]
  reg[i, "beta1"] <- sucrose_regression$coefficients[2]
  reg[i, "MSE"] <- sum(sucrose_regression$residuals^2) / sucrose_regression$df.residual
  reg[i, "Rsq"] <- summary(sucrose_regression)$r.squared
  reg[i, "pval"] <- summary(sucrose_regression)$coefficients[2,4]
}  

spread <- range(data.1[,15])
stripchart(data.1$Sucrose.dwb ~ data.1$Subsample, 
           vertical = TRUE, 
           main = "",
           xlab="",
           ylab="",
           type="n")
for(i in 1:nrow(reg)){
  abline(a=reg[i,1], b=(reg[i,2]))
}
abline(a=model3$coefficients[1], b=model3$coefficients[2], col="red", lwd=3)
text(x=2, y=20,"p-value: 0.00169 **")

## Regression on all samples for asparagine
reg <- matrix(nrow = length(genotypes), ncol = 5, dimnames = list(genotypes, c("beta0", "beta1", "MSE", "Rsq", "pval")))

for(i in 1:length(genotypes)){
  temp <- data.1[data.1$Genotype == genotypes[i],]
  asparagine_regression <- lm(temp[,4] ~ temp[,16])
  reg[i, "beta0"] <- asparagine_regression$coefficients[1]
  reg[i, "beta1"] <- asparagine_regression$coefficients[2]
  reg[i, "MSE"] <- sum(asparagine_regression$residuals^2) / asparagine_regression$df.residual
  reg[i, "Rsq"] <- summary(asparagine_regression)$r.squared
  reg[i, "pval"] <- summary(asparagine_regression)$coefficients[2,4]
}  

spread <- range(data.1[,4])
stripchart(data.1$Asparagine.dwb ~ data.1$Subsample, 
           vertical = TRUE, 
           main = "",
           xlab="",
           ylab="",
           type="n")
for(i in 1:nrow(reg)){
  abline(a=reg[i,1], b=(reg[i,2]))
}
abline(a=model4$coefficients[1], b=model4$coefficients[2], col="red", lwd=3)
text(x=2.25, y=3,"p-value: <2e-16 ***")

## Regression on all samples for starch dwb
reg <- matrix(nrow = length(genotypes), ncol = 5, dimnames = list(genotypes, c("beta0", "beta1", "MSE", "Rsq", "pval")))

for(i in 1:length(genotypes)){
  temp <- data.1[data.1$Genotype == genotypes[i],]
  starch_dwb_regression <- lm(temp[,13] ~ temp[,16])
  reg[i, "beta0"] <- starch_dwb_regression$coefficients[1]
  reg[i, "beta1"] <- starch_dwb_regression$coefficients[2]
  reg[i, "MSE"] <- sum(starch_dwb_regression$residuals^2) / starch_dwb_regression$df.residual
  reg[i, "Rsq"] <- summary(starch_dwb_regression)$r.squared
  reg[i, "pval"] <- summary(starch_dwb_regression)$coefficients[2,4]
}  

spread <- range(data.1[,13])
stripchart(data.1$Starch.dwb ~ data.1$Subsample, 
           vertical = TRUE, 
           main = "",
           xlab="",
           ylab="",
           type="n")
for(i in 1:nrow(reg)){
  abline(a=reg[i,1], b=(reg[i,2]))
}
abline(a=model6$coefficients[1], b=model6$coefficients[2], col="red", lwd=3)
text(x=2, y=76,"p-value: <2e-16 ***")

## Regression on all samples for crude protein
## Loading in the emmeans dataset
data.1 <- read.csv("data/all_emmeans_final.csv", header = T, stringsAsFactors = F)
#removing nej and ws
data.1 <- data.1[c(-361:-600),]

# deleting rows with missing data for crude protein
data.1 <- data.1[c(-90,-210,-330),]

## Setting variables as factors for genotype and subsample type
data.1[,1] <- as.factor(data.1[,1])
data.1[,2] <- as.factor(data.1[,2])
data.1[,2] <- factor(data.1[,2], levels = c("RK", "CK", "SK"))
str(data.1)

## Regression of all the data
model <- lm(Crude.Protein ~ Sample, data = data.1)
summary(model)

genotypes <- unique(data.1$Genotype)
reg <- matrix(nrow = length(genotypes), ncol = 5, dimnames = list(genotypes, c("beta0", "beta1", "MSE", "Rsq", "pval")))

for(i in 1:length(genotypes)){
  temp <- data.1[data.1$Genotype == genotypes[i],]
  protein_regression <- lm(temp[,5] ~ temp[,16])
  reg[i, "beta0"] <- protein_regression$coefficients[1]
  reg[i, "beta1"] <- protein_regression$coefficients[2]
  reg[i, "MSE"] <- sum(protein_regression$residuals^2) / protein_regression$df.residual
  reg[i, "Rsq"] <- summary(protein_regression)$r.squared
  reg[i, "pval"] <- summary(protein_regression)$coefficients[2,4]
}  

spread <- range(data.1[,5])
stripchart(data.1$Crude.Protein ~ data.1$Subsample, 
           vertical = TRUE, 
           main = "",
           xlab="",
           ylab="",
           type="n")
for(i in 1:nrow(reg)){
  abline(a=reg[i,1], b=(reg[i,2]))
}
abline(a=model$coefficients[1], b=model$coefficients[2], col="red", lwd=3)
text(x=1.5, y=14.25,"p-value: 0.288 ns")

### Figure 5 ###
## Clearing the global environment
rm(list=ls(all=TRUE))

## Loading in packages
library(ggplot2)
library(reshape2)
library(dplyr)
library(pheatmap)

# laoding in the data
data <- read.csv("data/all_emmeans_heatmap.csv", header=T, stringsAsFactors=F)
str(data)

# deleting columns that aren't needed
data <- data[,c(-1)]
str(data)
# compute correlation matrix
cormat <- round(cor(data, use="complete.obs"), 2)
melted_cormat <- melt(cormat)
# ggplot
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile()

#Get the upper triangle of the correlation matrix
get_upper_tri<-function(cormat){
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

## ggplot
ggheatmap <- ggplot(data = melted_cormat, aes(Var2, Var1, fill=value))+
  geom_tile(color="white") +
  scale_fill_gradient2(low="blue",high="red",mid="white",
                       midpoint = 0, limit=c(-1,1), space="Lab",
                       name="Pearson\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1)) +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label=value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1,0),
    legend.position = c(0.6,0.7),
    legend.direction = "horizontal") +
  guides(fill=guide_colorbar(barwidth = 7, barheight = 1,
                             title.position = "top", title.hjust = 0.5))