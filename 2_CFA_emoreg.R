###############################################################################
################## CFA for emotion regulation ER latent factor ################
###############################################################################

#=================================================================================
# CONTENTS
# 1. CFA of pre-specified emotion regulation model - primary measures
# 2. CFA of pre-specified emotion regulation model - primary + secondary measures
#==================================================================================

#----------------------------------------------------------------------------------
# Primary measures were: 
# amygdala-mPFC (dmPFC/vmPFC) resting state functional connectivity; 
# mPFC volume, basal forebrain volume, insula volume 
# lateral orbitofrontal cortical thickness

# Secondary measures were:
# LC-mPFC (vmPFC/dmPFC) and LC-amygdala resting state functional connectivity 
# DTI measures of LC-TEC bundle
# Negative reappraisal task score from Cam-CAN - scaled from negative to positive (i.e. higher scores=better reappraisal)

# All structural volumes are adjusted for (divided by) TIV
# All continuous data were scaled so that mean=5, SD=2
#----------------------------------------------------------------------------------

#-------------------------------------------------------
# 1. CFA pre-specified emoreg model - primary measures
#-------------------------------------------------------

emomodel <-'
emoreg=~BasCbr.FobBr_TIVavg + MedFroCbr_TIVavg + Ins_TIVavg+lateralorbitofrontal
+ dmPFC_amyg_avg + vmPFC_amyg_avg

BasCbr.FobBr_TIVavg~~MedFroCbr_TIVavg
BasCbr.FobBr_TIVavg~~Ins_TIVavg
MedFroCbr_TIVavg~~Ins_TIVavg
dmPFC_amyg_avg~~vmPFC_amyg_avg
'

model.fit <- cfa(emomodel, data=data, missing = "fiml", estimator="mlr") 
fitMeasures(model.fit, c("cfi","rmsea","srmr"))
MI<-modificationIndices(model.fit, sort=TRUE) 
MI
summary (model.fit, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)

#------------------------------------------------------
# 2. CFA pre-specified ER model with secondary measures 
#------------------------------------------------------

emomodel <-'
    emoreg=~ BasCbr.FobBr_TIVavg + MedFroCbr_TIVavg + Ins_TIVavg + 
       dmPFC_amyg_avg + vmPFC_amyg_avg +
       lateralorbitofrontal+
       meanMD_avg + meanRD_L+meanRD_R+meanFA_L+meanFA_R+
       dmPFC_LC + vmPFC_LC+LC_lamyg+LC_ramyg
      #+Neg_scale_reappraisal
         
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
  dmPFC_LC~~vmPFC_LC
  
  meanMD_avg~~meanRD_L
  meanMD_avg~~meanRD_R
  meanMD_avg~~meanFA_L
  meanMD_avg~~meanFA_R
  meanRD_L~~meanRD_R
  meanRD_L~~meanFA_L
  meanRD_L~~meanFA_R
  meanRD_R~~meanFA_L
  meanRD_R~~meanFA_R
  meanFA_L~~meanFA_R
    '
model.fit <- cfa(emomodel, data=data, missing = "fiml", estimator="mlr") 
fitMeasures(model.fit, c("cfi","rmsea","srmr"))
MI<-modificationIndices(model.fit, sort=TRUE)  #suggests accounting for covariance betwen dmPFC_LC amd vmPFC_LC, mean FA_L and meanFA_R
MI
summary (model.fit, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)
