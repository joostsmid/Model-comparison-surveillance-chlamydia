### R script to read in the key datasets for this project & to prepare them as R matrices

### Supposes existence of a subdirectory called 'data'

population.data <- read.csv("data/population_UK.csv")
m.15.19 <- population.data[,2]
m.20.24 <- population.data[,3]
m.25.34 <- population.data[,4]
m.35.44 <- population.data[,5]
f.15.19 <- population.data[,7]
f.20.24 <- population.data[,8]
f.25.34 <- population.data[,9]
f.35.44 <- population.data[,10]

nyears <- dim(population.data)[1]

if (dataset == "maxdata"){
  diagnoses.data <- read.csv("data/notificationS_UKmax.csv")
  tests.data <- read.csv("data/tests_UKmax.csv")

} else if (dataset == "mindata"){
  diagnoses.data <- read.csv("data/notificationS_UKmin.csv")
  tests.data <- read.csv("data/tests_UKmin.csv")
}

## notification data: per 100,000 people
notifications.m <- matrix(as.numeric(t(as.matrix(diagnoses.data[,2:5]*population.data[,2:5]/100000))), nrow = 4)
notifications.f <- matrix(as.numeric(t(as.matrix(diagnoses.data[,7:10]*population.data[,7:10]/100000))), nrow = 4)

## tests data: per 100 people
tested.m <- matrix(as.numeric(t(as.matrix(tests.data[,2:5]*population.data[,2:5]/100))), nrow = 4)
tested.f <- matrix(as.numeric(t(as.matrix(tests.data[,7:10]*population.data[,7:10]/100))), nrow = 4)

