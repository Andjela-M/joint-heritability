library(lavaan)
library(blavaan)
library(MASS)
library(tidySEM)
library(dplyr)
library(semTools)
library(ggplot2)

inputFileName = 'MRI_EEG.csv'

#read data
jointData <- read.csv(file=inputFileName)
jointData <- jointData[-2,] #remove mixed-sex twin pair
jointData$Status <- factor(jointData$Status, levels= c(1:2),labels= c("MZ", "DZ"))
dataset <- jointData[c("Status","Sigma1","Sigma2", "Thalamus1", "Thalamus2")]
  
# standardisation
EEG <- c(dataset$Sigma1, dataset$Sigma2)
EEG[EEG==-99] <- NA
EEG <- EEG %>% scale()
dataset$Sigma1 <- EEG[1:22]
dataset$Sigma2 <- EEG[23:44]
MRI <- c(dataset$Thalamus1, dataset$Thalamus2)
MRI[MRI==-99] <-  NA
MRI <- MRI %>% scale()
dataset$Thalamus1 <- MRI[1:22]
dataset$Thalamus2 <- MRI[23:44]
dataset <- na.omit(dataset)

colnames(dataset) <- c("zyg","P1", "P2", "P3", "P4")

#EEG
ace.modelEEG<-"
  A1=~ NA*P1 + c(a,a)*P1 
  A2=~ NA*P2 + c(a,a)*P2 
  C1 =~ NA*P1 + c(c,c)*P1
  C2 =~ NA*P2 + c(c,c)*P2
  # variances
  A1 ~~ 1*A1
  A2 ~~ 1*A2
  C1 ~~ 1*C1
  C2 ~~ 1*C2 
  P1~~c(e2,e2)*P1 
  P2~~c(e2,e2)*P2
  # covariances
  A1 ~~ c(1,.5)*A2 
  A1 ~~ 0*C1 + 0*C2 
  A2 ~~ 0*C1 + 0*C2 
  C1 ~~ c(1,1)*C2
  "
ace.fitEEGfreq <- lavaan(ace.modelEEG, data = dataset,group = "zyg")

# set priors, we use the defaults as they are reasonable with our 
# standardized data. 
# https://urldefense.com/v3/__https://ecmerkle.github.io/blavaan/articles/prior.html__;!!Dc8iu7o!ySweCoUdcPrwXp8ZLQ7CJZTynVlKv2GwScFCRRjOSM5PXY-YFBb98epxeCczv_ms0YW5M1i9ZY_GH1hHAibY$ 
priors <- dpriors(target = "stan")

## fit null model to calculate CFI, TLI, and NFI
null.model <- c(paste0("P", 1:2, " ~~ P", 1:2), paste0("P", 1:2, " ~ 1"))
fit0EEG <- blavaan(null.model, data = dataset, group = "zyg",
                   n.chains = 3, burnin = 1000, sample = 2000, seed = 20231009)

# fit blavaan
ace.fitEEG.blavaan <- blavaan::blavaan(ace.modelEEG,
                                       dp = priors,
                                       burnin = 1000,
                                       sample = 2000,
                                       seed = 20231009,
                                       data = dataset,
                                       group = "zyg",
                                       std.lv = TRUE)
summary(ace.fitEEG.blavaan)

# contributions with uncertainty
# extract stanfit object
stanfitEEG <- blavInspect(ace.fitEEG.blavaan,"mcobj")
saveRDS(stanfitEEG, file = paste0("Fit_channel", k, "_EEG.rds"))
drawsEEG <- posterior::as_draws_matrix(stanfitEEG)
draws_ACE <- drawsEEG[, 1:3] 
colnames(draws_ACE) <- c("A", "C", "E")
head(draws_ACE)
draws_ACE_prop <- draws_ACE %>% abs() %>% prop.table(., margin = 1) %>% 
  round(., 2)

#output
print(apply(draws_ACE_prop, 2, mean))
print(sum(apply(draws_ACE_prop, 2, mean)))

#model fit
fitIndices.fitEEG <- blavFitIndices(ace.fitEEG.blavaan, baseline.model = fit0EEG)
print(summary(fitIndices.fitEEG))

#plot
bayesplot::mcmc_dens(draws_ACE_prop) + xlim(0, 1) 
  
#MRI
ace.modelMRI<-"
  A3=~ NA*P3 + c(a,a)*P3 
  A4=~ NA*P4 + c(a,a)*P4 
  C3 =~ NA*P3 + c(c,c)*P3
  C4 =~ NA*P4 + c(c,c)*P4
  # variances
  A3 ~~ 1*A3
  A4 ~~ 1*A4
  C3 ~~ 1*C3
  C4 ~~ 1*C4 
  P3~~c(e2,e2)*P3 
  P4~~c(e2,e2)*P4
  # covariances
  A3 ~~ c(1,.5)*A4 
  A3 ~~ 0*C3 + 0*C4 
  A4 ~~ 0*C3 + 0*C4 
  C3 ~~ c(1,1)*C4
  "

ace.fitMRIfreq <- lavaan(ace.modelMRI, data = dataset,group = "zyg")

## fit null model to calculate CFI, TLI, and NFI
null.model <- c(paste0("P", 3:4, " ~~ P", 3:4), paste0("P", 3:4, " ~ 1"))
fit0MRI <- blavaan(null.model, data = dataset, group = "zyg",
                   n.chains = 3, burnin = 1000, sample = 2000, seed = 20231009)

# fit blavaan
ace.fitMRI.blavaan <- blavaan::blavaan(ace.modelMRI,
                                       dp = priors,
                                       burnin = 1000,
                                       sample = 2000,
                                       seed = 20231009,
                                       data = dataset,
                                       group = "zyg",
                                       std.lv = TRUE)
summary(ace.fitMRI.blavaan)
  
# contributions with uncertainty
# extract stanfit object
stanfitMRI <- blavInspect(ace.fitMRI.blavaan,"mcobj")
saveRDS(stanfitMRI, file = paste0("Fit_MRI.rds"))
drawsMRI <- posterior::as_draws_matrix(stanfitMRI)
draws_ACE <- drawsMRI[, 1:3] 
colnames(draws_ACE) <- c("A", "C", "E")
head(draws_ACE)
draws_ACE_prop <- draws_ACE %>% abs() %>% prop.table(., margin = 1) %>% 
  round(., 2)

#output
print(apply(draws_ACE_prop, 2, mean))
print(sum(apply(draws_ACE_prop, 2, mean)))

#model fit
fitIndices.fitMRI <- blavFitIndices(ace.fitMRI.blavaan, baseline.model = fit0MRI)
print(summary(fitIndices.fitMRI))

#plot
bayesplot::mcmc_dens(draws_ACE_prop) + xlim(0, 1) 

#EEG&MRI
ace.modelJoint<-"
  A1=~ NA*P1 + c(a,a)*P1 
  A2=~ NA*P2 + c(a,a)*P2 
  A3 =~ NA * P3 + c(a,a) * P3
  A4 =~ NA * P4 + c(a,a) * P4
  C1 =~ NA*P1 + c(c,c)*P1
  C2 =~ NA*P2 + c(c,c)*P2
  C3 =~ NA*P3 + c(c,c)*P3
  C4 =~ NA*P4 + c(c,c)*P4
  # variances
  A1 ~~ 1*A1
  A2 ~~ 1*A2
  A3 ~~ 1*A3
  A4 ~~ 1*A4
  C1 ~~ 1*C1
  C2 ~~ 1*C2 
  C3 ~~ 1*C3 
  C4 ~~ 1*C4 
  P1~~c(e2,e2)*P1 
  P2~~c(e2,e2)*P2
  P3~~c(e2,e2)*P3
  P4~~c(e2,e2)*P4
  # covariances
  A1 ~~ c(1,.5)*A2 
  A3 ~~ c(1,.5)*A4 
  A1 ~~ 0 * A3 + 0 * A4
  A2 ~~ 0 * A3 + 0 * A4
  A1 ~~ 0*C1 + 0*C2 + 0 *C3 + 0 * C4
  A2 ~~ 0*C1 + 0*C2 + 0 *C3 + 0 * C4
  A3 ~~ 0*C1 + 0*C2 + 0 *C3 + 0 * C4
  A4 ~~ 0*C1 + 0*C2 + 0 *C3 + 0 * C4
  C1 ~~ c(1,1)*C2
  C3 ~~ c(1,1)*C4
  C1 ~~ 0 * C3 + 0 * C4
  C2 ~~ 0 * C3 + 0 * C4
  "
  
## fit null model to calculate CFI, TLI, and NFI
null.model <- c(paste0("P", 1:4, " ~~ P", 1:4), paste0("P", 1:4, " ~ 1"))
fit0Joint <- blavaan(null.model, data = dataset, group = "zyg",
                     n.chains = 3, burnin = 1000, sample = 2000, seed = 20231009)

# fit blavaan
ace.fitJoint.blavaan <- blavaan::blavaan(ace.modelJoint,
                                         dp = priors,
                                         burnin = 1000,
                                         sample = 2000,
                                         seed = 20231009,
                                         data = dataset,
                                         group = "zyg",
                                         std.lv = TRUE)
summary(ace.fitJoint.blavaan)

# contributions with uncertainty
# extract stanfit object
stanfitJoint <- blavInspect(ace.fitJoint.blavaan,"mcobj")
saveRDS(stanfitJoint, file = paste0("Fit_channel", k, "_joint.rds"))
drawsJoint <- posterior::as_draws_matrix(stanfitJoint)
draws_ACE <- drawsJoint[, 1:3] 
colnames(draws_ACE) <- c("A", "C", "E")
head(draws_ACE)
draws_ACE_prop <- draws_ACE %>% abs() %>% prop.table(., margin = 1) %>% 
  round(., 2)

#output
print(apply(draws_ACE_prop, 2, mean))
print(sum(apply(draws_ACE_prop, 2, mean)))

#model fit
fitIndices.fitJoint <- blavFitIndices(ace.fitJoint.blavaan, baseline.model = fit0Joint)
print(summary(fitIndices.fitJoint))

#plot
bayesplot::mcmc_dens(draws_ACE_prop) + xlim(0, 1) 

#plot all distributions
d1 <- drawsEEG[, 1:3] %>% abs() %>% prop.table(., margin = 1) %>%
  round(., 2) %>% data.frame()
d2 <- drawsMRI[, 1:3] %>% abs() %>% prop.table(., margin = 1) %>%
  round(., 2)  %>% data.frame()
d3 <- drawsJoint[, 1:3] %>% abs() %>% prop.table(., margin = 1) %>%
  round(., 2)  %>% data.frame()
colnames(d1) <- c("A", "C", "E")
colnames(d2) <- c("A", "C", "E")
colnames(d3) <- c("A", "C", "E")
d1$Feature <- "EEG"
d2$Feature <- "MRI"
d3$Feature <- "EEG and MRI"
d <- rbind(d1, d2, d3)
long_data <- reshape2::melt(d, id.vars = "Feature")
df <- long_data
p <- ggplot(df, aes(x = value, fill = Feature)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ variable, scales = "free_x") +
  labs(title = paste0("Channel ", k),
       x = "Value",
       y = "Density") +
  xlim(0, 1) +
  theme_minimal() +
  theme(
    text = element_text(size = 17),
    axis.text = element_text(size = 13))
png(paste("Channel_", k, ".png", sep = ""), width=1000, height=800, res=120)
print(p)
dev.off()
