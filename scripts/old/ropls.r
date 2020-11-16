library(ropls)


rm(list=ls())


data <- read.table("proteinGroups_tryptic_reporterIntensities.csv", sep = ",", header = TRUE, row.names = 1)
dataMatrix <- data.matrix(data)
#dataMatrix[is.na(dataMatrix)] = 0


#dataMatrix <- dataMatrix[1:300,] #SUBSETTING
dataMatrix <- t(dataMatrix)
sampleMetadata <- read.table("proteinGroups_tryptic_sampleMetadata.csv", sep = ",", header = TRUE, row.names = 1)

# metadata
treatment <- sampleMetadata[, "treatment"]
states <- sampleMetadata[, "state"]
cell_line <- sampleMetadata[, "cell_line"]
replicate <- sampleMetadata[, "replicate"]

#sampleMetadata[]
#view(dataMatrix)
#sacurine.pca <- opls(dataMatrix)

data.pca <- opls(dataMatrix)
data.pca <- opls(dataMatrix, states)

plot(data.pca,
     typeVc = "x-score",
     parAsColFcVn = states)



####


# Select treatment
treatment_str <- "Nutlin"
dataMatrix_ <- data[,which(treatment %in% treatment_str)]
#dataMatrix_ <- t(dataMatrix_)

states_ <- states[which(treatment %in% treatment_str)]

data.pca_ <- opls(dataMatrix_)

plot(data.pca_,
     typeVc = "x-score",
     parAsColFcVn = states_)










