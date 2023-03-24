########################################
###### SEM - correlational models ###### 
########################################

#========================================================================================
# CONTENTS
# 1. Correlational model between resting HRV (ECG) and ER (primary + secondary measures)
# 1a. Multigroup analysis
# 2. Correlational model between HRV reactivity (PPG) and ER (primary+secondary measures)
# 2a. Measurement invariance testing for multigroup analysis
# 3. Exploratory correlational model between resting HRV and 4 latent emotionality factors 
# 3a. Multigroup analysis
# 4. Exploratory correlational model between HRV reactivity (PPG) and 4 latent emotionality factors
# 4a. Multigroup analysis

# age, where included, is shown in the models.
#========================================================================================

#------------------------------------------------------
# 1. resting HRV and primary/secondary ER + separate ERpsych
#------------------------------------------------------

model <-'

loghfrmssd_avg~HR # averaged HR/RMSSD HRV metric 

# neural ER measures
emoreg=~ BasCbr.FobBr_TIVavg + MedFroCbr_TIVavg + Ins_TIVavg + 
  dmPFC_amyg_avg + vmPFC_amyg_avg +
  lateralorbitofrontal+
  meanMD_avg + meanRD_L+meanRD_R+meanFA_L+meanFA_R+
  dmPFC_LC + vmPFC_LC+LC_lamyg+LC_ramyg

BasCbr.FobBr_TIVavg~~MedFroCbr_TIVavg
BasCbr.FobBr_TIVavg~~Ins_TIVavg
MedFroCbr_TIVavg~~Ins_TIVavg
dmPFC_amyg_avg~~vmPFC_amyg_avg
dmPFC_amyg_avg~~dmPFC_LC
dmPFC_amyg_avg~~vmPFC_LC
dmPFC_amyg_avg~~LC_lamyg
dmPFC_amyg_avg~~LC_ramyg
vmPFC_amyg_avg~~dmPFC_LC
vmPFC_amyg_avg~~vmPFC_LC
vmPFC_amyg_avg~~LC_lamyg
vmPFC_amyg_avg~~LC_ramyg
dmPFC_LC~~LC_lamyg
dmPFC_LC~~LC_ramyg
vmPFC_LC~~LC_lamyg
vmPFC_LC~~LC_ramyg
LC_lamyg~~LC_ramyg

meanMD_avg~~meanRD_L
meanMD_avg~~meanRD_R
meanMD_avg~~meanFA_L
meanMD_avg~~meanFA_R
meanRD_L~~meanRD_R
meanRD_L~~meanFA_L
meanRD_L~~meanFA_R
meanRD_R~~meanFA_L
meanRD_R~~meanFA_R
meanFA_L~~meanFA_L

meanFA_L~~meanFA_R
dmPFC_LC~~vmPFC_LC

#covariances
loghfrmssd_avg~~emoreg
#loghfrmssd_avg~~0*emoreg # for whole-group LRT to test significance of parameter

loghfrmssd_avg~~Neg_scale_reappraisal
#loghfrmssd_avg~~0*Neg_scale_reappraisal # for whole-group LRT to test significance of parameter

Neg_scale_reappraisal~~emoreg

#adjust for age
loghfrmssd_avg~age
emoreg~age
Neg_scale_reappraisal~age
HR~age
'
model.fit <- cfa(model, data=data, missing='FIML',estimator="mlr") #missing=ML
fitMeasures(model.fit, c("cfi","rmsea","srmr"))
summary (model.fit, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)

#compare freely estimated to constrained model
constrained <- cfa(model, data=data, missing='FIML',estimator="mlr") 
anova(model.fit,constrained) 

#-------------------------------------------------
# 1a. Multigroup analysis for [1]
#-------------------------------------------------
# Check age group invariance. Weak invariance achieved. 
config <- cfa(model, data=data, group="agegroup",missing = "fiml")
weak <- cfa(model, data=data, group="agegroup", group.equal=c("loadings"), missing = "fiml")
strong <- cfa(model, data=data, group="agegroup", group.equal=c("loadings", "intercepts"), missing = "fiml")
strict <- cfa(model, data=data, group="agegroup", group.equal=c("loadings", "intercepts", "residuals"), missing = "fiml")
structural <- cfa(model, data=data, group="agegroup", group.equal=c("loadings", "intercepts", "residuals", "means"), missing = "fiml")
#chi sq diff test
lavTestLRT(config, weak, strong, strict, structural) 

model <-'

loghfrmssd_avg~HR # averaged HR/RMSSD HRV metric 

# neural ER measures
emoreg=~ BasCbr.FobBr_TIVavg + MedFroCbr_TIVavg + Ins_TIVavg + 
  dmPFC_amyg_avg + vmPFC_amyg_avg +
  lateralorbitofrontal+
  meanMD_avg + meanRD_L+meanRD_R+meanFA_L+meanFA_R+
  dmPFC_LC + vmPFC_LC+LC_lamyg+LC_ramyg

BasCbr.FobBr_TIVavg~~MedFroCbr_TIVavg
BasCbr.FobBr_TIVavg~~Ins_TIVavg
MedFroCbr_TIVavg~~Ins_TIVavg
dmPFC_amyg_avg~~vmPFC_amyg_avg
dmPFC_amyg_avg~~dmPFC_LC
dmPFC_amyg_avg~~vmPFC_LC
dmPFC_amyg_avg~~LC_lamyg
dmPFC_amyg_avg~~LC_ramyg
vmPFC_amyg_avg~~dmPFC_LC
vmPFC_amyg_avg~~vmPFC_LC
vmPFC_amyg_avg~~LC_lamyg
vmPFC_amyg_avg~~LC_ramyg
dmPFC_LC~~LC_lamyg
dmPFC_LC~~LC_ramyg
vmPFC_LC~~LC_lamyg
vmPFC_LC~~LC_ramyg
LC_lamyg~~LC_ramyg

meanMD_avg~~meanRD_L
meanMD_avg~~meanRD_R
meanMD_avg~~meanFA_L
meanMD_avg~~meanFA_R
meanRD_L~~meanRD_R
meanRD_L~~meanFA_L
meanRD_L~~meanFA_R
meanRD_R~~meanFA_L
meanRD_R~~meanFA_R
meanFA_L~~meanFA_L

meanFA_L~~meanFA_R
dmPFC_LC~~vmPFC_LC

#covariances
loghfrmssd_avg~~emoreg
loghfrmssd_avg~~Neg_scale_reappraisal
Neg_scale_reappraisal~~emoreg

#multigroup LRT comparisons for covariances of interest
#loghfrmssd_avg~~c(young_effect,old_effect)*emoreg # LRT test. freely estimated model
#loghfrmssd_avg~~c(effect, effect)*emoreg  #constrained to be equal model
#loghfrmssd_avg~~c(0,NA)*emoreg # LRT test for if path is significant in either group

#loghfrmssd_avg~~c(young_effect,old_effect)*Neg_scale_reappraisal   # LRT test. freely estimated model
#loghfrmssd_avg~~c(effect,effect)*Neg_scale_reappraisal  #constrained to be equal model
#loghfrmssd_avg~~c(0,NA)*Neg_scale_reappraisal # LRT test for if path is significant in either group

#adjust for age
loghfrmssd_avg~age
emoreg~age
Neg_scale_reappraisal~age
HR~age
'

#age group analysis younger vs older
fit.age <- cfa(model, data=data, group="agegroup", group.equal=c("loadings"), missing = "fiml",estimator="mlr") 
summary(fit.age, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)
fitMeasures(fit.age, c("cfi","rmsea","srmr"))

#LRT comparing freely estimated to constrained models 
constrained <- cfa(model, data=data, group="agegroup", group.equal = c("loadings"), missing = "fiml", estimator="mlr") #missing=ML
anova(fit.age,constrained) 

#========================================================================
# 2. HRV reactivity from PPG DATA and 'neural' ER and Neg reappraisal
#=========================================================================
model<-'
#PPG latent change score model
loghfrmssdmov_PPGavg ~ 1*loghfrmssd_PPGavg     # Fixed regression of loghfrmssdmov_PPGavg on loghfrmssd_PPGavg
dhrv1 =~ 1*loghfrmssdmov_PPGavg     # Fixed regression of dhrv1 on loghfrmssdmov_PPGavg
loghfrmssdmov_PPGavg ~ 0*1          # This line constrains the intercept of loghfrmssdmov_PPGavg to 0
loghfrmssdmov_PPGavg ~~ 0*loghfrmssdmov_PPGavg    # This fixes the variance of the loghfrmssdmov_PPGavg to 0 

dhrv1 ~ 1             # This estimates the intercept of the change scores 
loghfrmssd_PPGavg ~  1           # This estimates the intercept of loghfrmssd_PPGavg 
dhrv1 ~~  dhrv1       # This estimates the variance of the change scores 
loghfrmssd_PPGavg ~~   loghfrmssd_PPGavg    # This estimates the variance of loghfrmssd_PPGavg 
dhrv1~loghfrmssd_PPGavg          # This estimates the self-feedback parameter

loghfrmssd_PPGavg~HR_rest
loghfrmssdmov_PPGavg~HR_mov
HR_mov~~HR_rest

#ER neural primary and secondary measures
emoreg=~ BasCbr.FobBr_TIVavg + MedFroCbr_TIVavg + Ins_TIVavg + 
  dmPFC_amyg_avg + vmPFC_amyg_avg +
  lateralorbitofrontal+
  meanMD_avg + meanRD_L+meanRD_R+meanFA_L+meanFA_R+
  dmPFC_LC + vmPFC_LC+LC_lamyg+LC_ramyg

BasCbr.FobBr_TIVavg~~MedFroCbr_TIVavg
BasCbr.FobBr_TIVavg~~Ins_TIVavg
MedFroCbr_TIVavg~~Ins_TIVavg
dmPFC_amyg_avg~~vmPFC_amyg_avg
dmPFC_amyg_avg~~dmPFC_LC
dmPFC_amyg_avg~~vmPFC_LC
dmPFC_amyg_avg~~LC_lamyg
dmPFC_amyg_avg~~LC_ramyg
vmPFC_amyg_avg~~dmPFC_LC
vmPFC_amyg_avg~~vmPFC_LC
vmPFC_amyg_avg~~LC_lamyg
vmPFC_amyg_avg~~LC_ramyg
dmPFC_LC~~LC_lamyg
vmPFC_LC~~LC_lamyg
dmPFC_LC~~LC_ramyg
vmPFC_LC~~LC_ramyg
LC_lamyg~~LC_ramyg

meanMD_avg~~meanRD_R
meanMD_avg~~meanRD_L
meanMD_avg~~meanFA_L
meanMD_avg~~meanFA_R
meanRD_L~~meanRD_R
meanRD_L~~meanFA_L
meanRD_L~~meanFA_R
meanRD_R~~meanFA_L
meanRD_R~~meanFA_R
meanFA_L~~meanFA_L

meanFA_L~~meanFA_R
dmPFC_LC~~vmPFC_LC

#covariances
emoreg~~dhrv1
Neg_scale_reappraisal~~dhrv1
Neg_scale_reappraisal~~emoreg

#dhrv1~~0*emoreg # constrained model for LRT
#dhrv1~~0*Neg_scale_reappraisal # constrained model for LRT

'

fit.model<- cfa(model, data=data, estimator='mlr',fixed.x=FALSE,missing='fiml')
summary(fit.model, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)
fitMeasures(fit.model, c("cfi","rmsea","srmr"))
constrained <- cfa(model, data=data, missing='FIML',estimator="mlr",fixed.x=FALSE) # HRV-ER path constrained to = 0
anova(fit.model,constrained) 
onyx(fit.model)

#--------------------------------------------------------------------
# 2a. Measurement invariance testsing for multigroup analysis for [2]
#--------------------------------------------------------------------
# configural invariance achieved. Metric invariance not achieved. Therefore multigroup analysis not conducted.

config <- cfa(model, data=data, group="agegroup", missing = "fiml")
weak <- cfa(model, data=data, group="agegroup", group.equal=c("loadings"), missing = "fiml")
strong <- cfa(model, data=data, group="agegroup", group.equal=c("loadings", "intercepts"), missing = "fiml")
strict <- cfa(model, data=data, group="agegroup", group.equal=c("loadings", "intercepts", "residuals"), missing = "fiml")
structural <- cfa(model, data=data, group="agegroup", group.equal=c("loadings", "intercepts", "residuals", "means"), missing = "fiml")
#chi sq diff test
lavTestLRT(config, weak, strong) 


#----------------------------------------------------
# 3. Resting HRV and 4 latent emotionality factors 
#----------------------------------------------------
model<-
  '
loghfrmssd_avg~HR

#4 latent ER factors
basneg_aff=~NeuW_mean_neg+PosW_mean_neg+NegW_mean_neg+PosW_mean_pos
posreac=~PosW_mean_pos+NeuW_mean_pos+NegW_mean_pos
posreg=~NegW_mean_pos+NegR_mean_pos
negreac=~NegW_mean_neg+NegR_mean_neg

PosW_mean_pos~~PosW_mean_neg
NegW_mean_pos~~NegR_mean_neg

# covariances
basneg_aff~~loghfrmssd_avg
posreac~~loghfrmssd_avg
posreg~~loghfrmssd_avg
#posreg~~0*loghfrmssd_avg #LRT test constrained path
negreac~~loghfrmssd_avg

# age adjustment
basneg_aff~age
posreac~ age
posreg~ age
negreac~age
loghfrmssd_avg~age
HR~age

#LC covariate
#basneg_aff~meRLC
#posreac~meRLC
#posreg~meRLC
#negreac~meRLC
#loghfrmssd_avg~meRLC
#meRLC~age
'
model.fit <- cfa(model, data=data, missing='FIML',estimator="mlr") #missing=ML
summary (model.fit, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)
fitMeasures(model.fit, c("cfi","rmsea","srmr"))
MI<-modificationIndices(model.fit, sort=TRUE)
MI
constrained <- cfa(model, data=data, missing='FIML',estimator="mlr") 
anova(model.fit,constrained) 

#----------------------------------------------------
# 3a. Multigroup analysis 
#----------------------------------------------------
# Check age group invariance. Weak invariance achieved. 
config <- cfa(model, data=data, group="agegroup",missing = "fiml")
weak <- cfa(model, data=data, group="agegroup", group.equal=c("loadings"), missing = "fiml")
strong <- cfa(model, data=data, group="agegroup", group.equal=c("loadings", "intercepts"), missing = "fiml")
strict <- cfa(model, data=data, group="agegroup", group.equal=c("loadings", "intercepts", "residuals"), missing = "fiml")
structural <- cfa(model, data=data, group="agegroup", group.equal=c("loadings", "intercepts", "residuals", "means"), missing = "fiml")
#chi sq diff test
lavTestLRT(config, weak, strong, strict, structural) 

model<-
  ' 
  loghfrmssd_avg~HR

#4 latent ER factors
basneg_aff=~NeuW_mean_neg+PosW_mean_neg+NegW_mean_neg+PosW_mean_pos
posreac=~PosW_mean_pos+NeuW_mean_pos+NegW_mean_pos
posreg=~NegW_mean_pos+NegR_mean_pos
negreac=~NegW_mean_neg+NegR_mean_neg

PosW_mean_pos~~PosW_mean_neg
NegW_mean_pos~~NegR_mean_neg

# covariances
basneg_aff~~loghfrmssd_avg
posreac~~loghfrmssd_avg
posreg~~loghfrmssd_avg
negreac~~loghfrmssd_avg

#multigroup comparisons
#posreg~~c(effect, effect)*loghfrmssd_avg #compares groups
#posreg~~c(0,NA)*loghfrmssd_avg # tests significance of path in either group

# age adjustment
basneg_aff~age
posreac~ age
posreg~ age
negreac~age
loghfrmssd_avg~age
HR~age

#LC covariate
#basneg_aff~meRLC
#posreac~meRLC
#posreg~meRLC
#negreac~meRLC
#loghfrmssd_avg~meRLC
#meRLC~age
'

#age group analysis younger vs older
fit.age <- cfa(model, data=data, group="agegroup", group.equal=c("loadings"), missing = "fiml",estimator="mlr") 
summary(fit.age, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)
fitMeasures(fit.age, c("cfi","rmsea","srmr"))

#LRT comparing freely estimated to constrained models 
constrained <- cfa(model, data=data, group="agegroup", group.equal = c("loadings"), missing = "fiml", estimator="mlr") #missing=ML
anova(fit.age,constrained) 

#----------------------------------------------------------
# 4. HRV reactivity (PPG) and 4 emotionality factors
#----------------------------------------------------------
model<-'
#PPG latent change score model
loghfrmssdmov_PPGavg ~ 1*loghfrmssd_PPGavg     # Fixed regression of loghfrmssdmov_PPGavg on loghfrmssd_PPGavg
dhrv1 =~ 1*loghfrmssdmov_PPGavg     # Fixed regression of dhrv1 on loghfrmssdmov_PPGavg
loghfrmssdmov_PPGavg ~ 0*1          # This line constrains the intercept of loghfrmssdmov_PPGavg to 0
loghfrmssdmov_PPGavg ~~ 0*loghfrmssdmov_PPGavg    # This fixes the variance of the loghfrmssdmov_PPGavg to 0 

dhrv1 ~ 1             # This estimates the intercept of the change scores 
loghfrmssd_PPGavg ~  1           # This estimates the intercept of loghfrmssd_PPGavg 
dhrv1 ~~  dhrv1       # This estimates the variance of the change scores 
loghfrmssd_PPGavg ~~   loghfrmssd_PPGavg    # This estimates the variance of loghfrmssd_PPGavg 
dhrv1~loghfrmssd_PPGavg          # This estimates the self-feedback parameter

loghfrmssd_PPGavg~HR_rest
loghfrmssdmov_PPGavg~HR_mov
HR_rest~~HR_mov

#4 latent ER factors
basneg_aff=~NeuW_mean_neg+PosW_mean_neg+NegW_mean_neg+PosW_mean_pos
posreac=~PosW_mean_pos+NeuW_mean_pos+NegW_mean_pos
posreg=~NegW_mean_pos+NegR_mean_pos
negreac=~NegW_mean_neg+NegR_mean_neg

PosW_mean_pos~~PosW_mean_neg
NegW_mean_pos~~NegR_mean_neg

basneg_aff~~dhrv1
posreac~~dhrv1
posreg~~dhrv1
negreac~~dhrv1
'

fit.model<- cfa(model, data=data, estimator='mlr',fixed.x=FALSE,missing='fiml')
summary(fit.model, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)
fitMeasures(fit.model, c("cfi","rmsea","srmr"))
constrained <- cfa(model, data=data, missing='FIML',estimator="mlr",fixed.x=FALSE) # HRV-ER path constrained to = 0
anova(fit.model,constrained) 

#----------------------------------------------------
# 4a. Multigroup analysis 
#----------------------------------------------------
# Check age group invariance. Weak metric invariance achieved. 
config <- cfa(model, data=data, group="agegroup",missing = "fiml")
weak <- cfa(model, data=data, group="agegroup", group.equal=c("loadings"), missing = "fiml")
strong <- cfa(model, data=data, group="agegroup", group.equal=c("loadings", "intercepts"), missing = "fiml")
strict <- cfa(model, data=data, group="agegroup", group.equal=c("loadings", "intercepts", "residuals"), missing = "fiml")
structural <- cfa(model, data=data, group="agegroup", group.equal=c("loadings", "intercepts", "residuals", "means"), missing = "fiml")
#chi sq diff test
lavTestLRT(config, weak, strong, strict, structural) 

fit.age <- cfa(model, data=data, group="agegroup", group.equal = c("loadings"), missing = "fiml", estimator="mlr") #missing=ML
summary(fit.age, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)
fitMeasures(fit.age, c("cfi","rmsea","srmr"))


