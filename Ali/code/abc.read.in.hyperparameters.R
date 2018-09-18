### R script to read in our hyperparameters for prior specification

### Supposes existence of a subdirectory called 'data'

### Required package
# library(xlsx) # (available from online R archive: CRAN)

# Use a different library for reading in excel files so rJava is 
# not required
# load.library(readxl) 

## Read in hyperparameters

# my file: variable definitions as indicated by row names in .xlsx file (& by usage below, of course!)
# hyperparameters <- read.xlsx("data/priors.xlsx",sheetIndex=1,colIndex=c(2:4),rowIndex=c(2:18),header=T) 

# hyperparameters <- read_excel("data/priors.xlsx", skip = 1) 
hyperparameters <- read.csv("data/priors.csv")
hyperparameters <- hyperparameters[,2:ncol(hyperparameters)]

prior.mvn.inf <- c(as.numeric(hyperparameters[10,1]),as.numeric(hyperparameters[11,1]),as.numeric(hyperparameters[12,1]),as.numeric(hyperparameters[13,1]))
prior.mvn.screening <- c(as.numeric(hyperparameters[14,1]),as.numeric(hyperparameters[15,1]),as.numeric(hyperparameters[16,1]),as.numeric(hyperparameters[17,1]))


