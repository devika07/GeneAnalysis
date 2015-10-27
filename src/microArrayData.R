###################################################################################################################################
# PROJECT FOR GENE EXPRESSION ANALYSIS
# TOPIC 1
###################################################################################################################################

#TASK 1

#Reads a file in table format and creates a data frame from it
#chosen separator of the file: the tab '\t'
#emit NA strings (NA or null values)
dataset <- read.table("GDS4879-ORIGINAL-withoutComments.soft", header=TRUE, sep="\t", strip.white=TRUE, na.strings=c("NA","null"))

# list rows of data that have missing values 
dataset[!complete.cases(dataset[1:10000,]),]

# create new dataset without missing data 
newdataset <- na.omit(dataset)

#choose a subset from the dataset 1-10000 lines:
newdataset <- newdataset[1:10000,]

#Save the second column into a variable gene.names
gene.names <- newdataset['IDENTIFIER']

#Select desired genes:
my.genes <- c("GSM1085677","GSM1085681","GSM1085685", "GSM1085689", "GSM1085695", "GSM1085698",
            "GSM1085673", "GSM1085679", "GSM1085694", "GSM1085696", "GSM1085699",
            "GSM1085701",
            "GSM1085666", "GSM1085668", "GSM1085670", "GSM1085671", "GSM1085674", "GSM1085678",
            "GSM1085665", "GSM1085667", "GSM1085669", "GSM1085672", "GSM1085675",
            "GSM1085676" )

#Final dataset consists of a subset of 24 samples:
final.dataset <- newdataset[my.genes]

###################################################################################################################################

# TASK 2:

# Boxplot of the dataset: 
boxplot(final.dataset,            #dataframe
        las=2,                    #vertical names
        par(mar = c(7, 2, 1, 1)), #increase size of the plot in the window
        col = c("red", "red", "red", "red", "red", "red",
                "palevioletred","palevioletred","palevioletred","palevioletred","palevioletred","palevioletred",
                "black", "black","black","black","black","black",
                "grey", "grey","grey","grey","grey","grey") #color on the box plot
        )

####################################################################################################################################

# TASK 3:

#Convert data from a frame format to a matrix format

gnames <- gene.names[,1]                # assign labels in column 1 to "rnames"
mat_data <- data.matrix(final.dataset[,1:ncol(final.dataset)])  # transform columns into a matrix
rownames(mat_data) <- gnames            # assign row names 

#
# Construct a heatmap of all samples
#
png("heatmap1.png", units = "px", width=8*640, height=8*640, res=300, pointsize=4)
heatmap1 <- heatmap(mat_data,          #matrix data
                    cexCol=0.7,        #size of col names
                    cexRow=0.2,        #size of row names
                    margins=c(3,2)
                   # scale="none"
                  )
dev.off()

#
# Construct a heatmap of all samples using the function heatmap.2 from the library gplots. 
# Use scale = 'none' and trace = 'none'. 
#

# Install and load if require the package for gplots
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}

# Install and load if required the RColorBrewer
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

#call heatmap2
png("heatmap2.png", units = "px", width=3840, height=2160, res=1024)
heatmap.2(mat_data,
          col=heat.colors(256), 
          scale="none", 
          key=T, 
          keysize=1.5,
          density.info="none", 
          trace="none",
          cexCol=0.9, 
          labRow=NA,
          margins =c(4,2)       # widens margins around plot
          )
dev.off()


####################################################################################################################################

# TASK 4
#make the dataset transponse
transposed.matrix <- data.matrix(final.dataset[0:nrow(final.dataset), 0:ncol(final.dataset)])
transposed.matrix <- t(transposed.matrix)


#make the pca analysis
pca <- prcomp(transposed.matrix, center = TRUE, scale = TRUE)
summary(pca)

#make the plot
plot(pca$x[, 1:2],  col = c("green", "blue", "red", "black"))

#export plot to high resolution .png image
png("pca_analysis_task4.png", units = "px", width=3840, height=2160, res=300)
plot(pca$x[, 1:2],  col = c("green", "blue", "red", "black"))
dev.off()

####################################################################################################################################

#TASK 5
#create empty vector
vec <- vector()
#calculate p-value for all rows of datase
for(i in 1:nrow(final.dataset)){
  males <- final.dataset[i, 1:12]
  females <- final.dataset[i, 13:24]
  
  ttest <- t.test(males,females)
  vec[i] <- ttest$p.value
}

#assign gene names to results
data = data.matrix(vec)
row.names(data) <- gnames

#sort the results
data.sort<-apply(data,2,sort)

#write results tou .csv file
write.table(data.sort[0:100, 0:ncol(data.sort)], file = "task5_pvalues.csv", sep = ",",
            col.names = FALSE, row.names = TRUE,
            qmethod = "double")

####################################################################################################################################

#TASK 6
#create empty vector
vec <- vector()
#calculate p-value for all rows of datase
for(i in 1:nrow(final.dataset)){
  alcoholics <- final.dataset[i, 1:6] + final.dataset[i, 13:18]
  non_alcoholics <- final.dataset[i, 7:12] + final.dataset[i, 19:24]
  
  ttest <- t.test(alcoholics, non_alcoholics)
  vec[i] <- ttest$p.value
}

#assign gene names to results
data = data.matrix(vec)
row.names(data) <- gnames

#sort the results
data.sort<-apply(data,2,sort)

#write results tou .csv file
write.table(data.sort[0:100, 0:ncol(data.sort)], file = "task6_pvalues.csv", sep = ",",
            col.names = FALSE, row.names = TRUE,
            qmethod = "double")

####################################################################################################################################
