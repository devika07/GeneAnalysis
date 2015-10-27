###################################################################################################################################
# PROJECT FOR GENE EXPRESSION ANALYSIS
# TOPIC 2 - from RNA Seq Dataset
###################################################################################################################################

dataset <- read.table(file.choose(),header=TRUE,stringsAsFactor=FALSE)
exclude <-which(dataset$significant!="yes" | dataset$value_1 > dataset$value_2)
dataset <- dataset[-exclude,]
