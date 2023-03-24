##################################################################
################## CFA for HRV latent factor ####################
##################################################################
library(lavaan)

#continuous data was scaled to have mean=5 and SD=2

#====================================================================================
                #CONTENTS#
# 1. CFA for theory-driven HRV factor (pre-registered)
# 2. CFA for exploratory bifactor HRV model
# 3  CFA for exploratory 2nd order HRV model 
# 4. CFA for exploratory resting HRV factor comprising averaged HF/RMSSD values (ECG)
# 5. CFA for exploratory HRV reactivity factor using averaged HF/RMSSD values in a latent score change model (PPG)

# 4. and 5. were subsequently used in the correlational models. 
#====================================================================================

#-----------------------------------
# 1. CFA for theory-driven HRV factor 
#-----------------------------------

HRVmodel <-'
HRV=~ logRMSSD+loghf+logRMSSD_diff+loghf_diff
logRMSSD~HR
loghf~HR
logRMSSD_diff~HRmean_diff
loghf_diff~HRmean_diff
HR~~HRmean_diff

logRMSSD~~loghf
loghf_diff~~logRMSSD_diff
'
model.fit <- cfa(HRVmodel, data=data, missing = "fiml", estimator="mlr") 
fitMeasures(model.fit, c("cfi","rmsea","srmr"))
MI<-modificationIndices(model.fit, sort=TRUE)
MI
summary (model.fit, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)

#-----------------------------
#2. Try bifactor model 'methods'
#-----------------------------
# Non-optimal fit; negative variances 

HRVmodel <-'
HRV=~ logRMSSD+loghf+logRMSSD_diff+loghf_diff
ECG=~logRMSSD+loghf
PPG=~logRMSSD_diff+loghf_diff
ECG~HR
PPG~HRmean_diff
ECG~~0*HRV #methods factors orthogonal to latent HRV factor
PPG~~0*HRV #methods factors orthogonal to latent HRV factor
ECG~~PPG
HR~~HRmean_diff
'
#-----------------------------
#3. Try second order model
#-----------------------------
#similar errors to above; no correlation between PPG and ECG first order factors

HRVmodel <-'
HRV=~ECG+PPG
ECG=~logRMSSD+loghf
PPG=~logRMSSD_diff+loghf_diff
ECG~HR
PPG~HRmean_diff
ECG~~PPG
HR~~HRmean_diff
'

#--------------------------------------------------------------------------------
# 4. A latent HRV factor comprising averaged HF/RMSSD values for resting ECG 
#--------------------------------------------------------------------------------

HRVmodel <-'
HRV=~loghfrmssd_avg+loghfrmssd_PPGavg+loghfrmssdmov_PPGavg
HRV~HR
loghfrmssd_PPGavg~~loghfrmssdmov_PPGavg
'
model.fit <- cfa(HRVmodel, data=data, missing = "fiml", estimator="mlr")
fitMeasures(model.fit, c("cfi","rmsea","srmr"))
summary (model.fit, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)
MI<-modificationIndices(model.fit, sort=TRUE)
MI

#------------------------------------------------------------------------
# 5. LATENT CHANGE SCORE MODEL WITH MOVIE/REST PULSE OX DATA - avg HF/RMSSD PPG only
# t1 (loghfrmssd_PPGavg) and t2 (loghfrmssdmov_PPGavg)
#------------------------------------------------------------------------

model<-'

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
'
model.fit <- cfa(model, data=data, estimator='mlr',fixed.x=FALSE,missing='fiml')
summary(model.fit, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)
fitMeasures(model.fit, c("cfi","rmsea","srmr"))
MI<-modificationIndices(model.fit, sort=TRUE)
MI

