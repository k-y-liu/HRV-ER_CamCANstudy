# HRV-ER_CamCANstudy
# SEM analyses

This repository contains the R code used to analyze data, using structural equation modelling, obtained from the Cambridge Centre for Ageing and Neuroscience (Cam-CAN) data repository (http://www.mrc-cbu.cam.ac.uk/datasets/camcan/). Mean, normalized locus coeruleus signal intensity (LC CR) values from participants' MT-weighted MRI images were obtained as described in this paper (https://www.sciencedirect.com/science/article/pii/S0197458018303786?via%3Dihub) and can be made available upon reasonable request. 

## Running the code

You will need to install [R](http://cran.rstudio.com/), [Rstudio](http://www.rstudio.com/ide/download/) and Lavaan. To install Lavaan, open Rstudio and run the commands:


```r
install.packages("lavaan", dependencies=TRUE)
library(lavaan)
```

Next, you need to import data by changing the working directory to the data you would like to analyze. 
