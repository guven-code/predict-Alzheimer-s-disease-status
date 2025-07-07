#### Alzhemier's dementia classifying (AI/ML) model prediction 
#modeling of biomarkers for the most correlated genes
#based on SET1 pTau/Abeta42 levels by gene expression
#and SET1 is partitioned into train(0.8) and test(0.2)
#later tested on SET2
##AUTHORs: EGUVEN and RMELLER
##UPDATED: 05/21/2025
## 2 different datasets REF https://datascience.stackexchange.com/questions/92996/how-to-improve-machine-learning-model-using-2-datasets
#### REFS https://topepo.github.io/caret/index.html
## REFS #https://www.machinelearningplus.com/machine-learning/caret-package/

#####LOAD THE PACKAGES
library(colorRamp2)
library(RColorBrewer)
library(statmod)
library(sva)
library(caret)
library(MatchIt)
library(tidyr)
library(ggplot2)
library(limma)
library(edgeR)
library(MLmetrics)
library(caretEnsemble)
library(dplyr)
library(ggfortify)
library(MLeval)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(enrichplot)
library(gridExtra)
library(corrplot)
library(ComplexHeatmap)
library(grid)
library(patchwork)  
library(glmnet)
library(boot)
library(ggpubr)
library(pROC)
library(Metrics)  
library(verification)
library(VennDiagram)


rm(list = ls())
gc() #free up memory and report memory usage

################################################################################
##### STEP 1:  Import data and clean up tables.           ######################
################################################################################
setwd("/Users/rmeller/Desktop/DESKTOP/ALZHEIMERS/FINAL_PRED")
getwd()
dir.create("Plots")
dir.create("RESULTS")


df_0 <- read.csv("gene_count_matrix.csv", header = TRUE)
##Load Datasheet (Patient Phenotype data from Emory-ADRC)
pheno<- read.table("Copy of AD_Patient_data_Set1_Set2_07092024_clean.txt", sep = "\t", header = TRUE)

### Use this to change the cutoff used to call AD status.  
## Schindler et al 2018 used a cutoff of 0.0198 vs 0.025 from Emory ADRC
## This will recode pTau/Ab to >0.0198
#pheno$CAT <- NULL
#pheno$CAT <- ifelse(pheno$"ptau.abeta" > 0.0198, "YES", "NO")

#head(df_0); 
dim (df_0) # 62031 x 169
#head(pheno); 
dim (pheno) # 168 x 17

## Remove samples which are repeats.  
pheno <- pheno[!pheno$REPEAT=="y",]
dim(pheno)# 164 samples
# remove any where the CAT=NA
pheno <- pheno[!is.na(pheno$ABETA),]
dim(pheno)# 154 samples
# change GUIDD tofrom 1.XXX to 1_XXX
pheno$GUID <- gsub("[.]","_", pheno$GUID)

## Modify and filter counts to creat correlation cpm data

#REMOVE STRINGTIE IDs sometimes separate function throws this error:
#Error: vector memory limit of 100.0 Gb reached, see mem.maxVSize()
exp_counts <- separate(df_0, col=gene_id, into= c("stringtie_id", "gene_id"),
                       extra = "merge", sep ="\\|", fill = "left")

# remove col stringtie_id
exp_counts$stringtie_id <- NULL
head(exp_counts,1); dim(exp_counts)


###MATCH THE PHENO_DATA SAMPLES WITH EXP_COUNTS
# change GUIDD tofrom 1.XXX to 1_XXX 
# IDs batch 1 start with an "X" so sub that out.
# IDs batch 2 start with nAD_ so sub that out
colnames(exp_counts) <- gsub("[.]","_", colnames(exp_counts))
colnames(exp_counts) <- gsub("[-]","_", colnames(exp_counts))
colnames(exp_counts) <- gsub("X","", colnames(exp_counts))
colnames(exp_counts) <- gsub("nAD_","", colnames(exp_counts))

# we have one repeated sample but with diff sex,  GUID 1.65166, so remove
pheno <- pheno[!pheno$GUID=="1_65166",]
dim(pheno)# 152 samples

# Now filter the IDS
exp_counts2 <- exp_counts[,names(exp_counts) %in% pheno$GUID]

pheno_data<-pheno[pheno$GUID %in% names(exp_counts2),]
dim(exp_counts2);# 62031 x 148
dim(pheno_data);# 148 x 17 

# Check and filter the training data
# now add back the Gene IDs
rownames(exp_counts2) <- make.unique(exp_counts$gene_id)
dim(exp_counts2);# 62031 x 148
dim(pheno_data);# 148 x 18 


## NOW let's look at the data we have. 
plot_pTau<-ggplot(pheno, aes(x = AGE, y = (pTAU.ABETA), ymin =0.003, ymax=0.3)) +
  geom_point(aes(color = SEX),size=2.5) +
  facet_wrap(~ SEX) +
  labs(title = "pTau/Abeta42 ratios vs. Age by Sex",
       x = "Age",
       y = "pTau/Abeta42") +
  scale_y_log10(breaks=c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1), labels=c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1)) +
  geom_hline(yintercept=(0.025), linetype="solid", color="black", linewidth=1) +
  theme_minimal()+theme(plot.title = element_text(face = "bold",size=14),
                        axis.title = element_text(face = "bold", size=14),
                        axis.text = element_text(size=14),
                        strip.text = element_text(face = "bold",size=14),
                        legend.title = element_text(face = "bold",size=14),
                        legend.text = element_text(size=14))

tiff(file = "Plots/01.PTAUabeta_ALL.tiff", width = 10, height = 8, units = "in", res = 300)
print(plot_pTau)
dev.off()

pheno$pTAU.ABETA

## Modify and filter counts to create correlation cpm data

## identify genes with low counts and store names
exp_counts3 <- exp_counts2[(rowSums(exp_counts2) > 150),]
dim(exp_counts3)# 20012 x 148
filt_genes <- rownames(exp_counts3)


# Step 4: Convert to CPM (Counts Per Million) for each dataset
cpm_genes <- apply(exp_counts2, 2, function(x) (x / sum(x)) * 1e6)
cpm_genes[1:10,1:10]
cat("CPM data dimensions:", dim(cpm_genes), "\n")

venn.diagram(
  x = list(colnames(cpm_genes), filt_genes),
  category.names = c("Set 1" , "Set 2 "),
  filename = 'Plots/#14_venn_diagramm.png',
  output=TRUE
)

## Now refine lists
cpm_genes <- cpm_genes[rownames(cpm_genes) %in% filt_genes,]
cat("CPM data dimensions:", dim(cpm_genes), "\n")
#cpm_genes <-t(cpm_genes)

## Split data into SET1 and SET2
SET1_pheno<-pheno_data[pheno_data$SET=="1",]
SET1_exp_counts<-cpm_genes[,pheno_data$SET=="1"]
dim(SET1_exp_counts)# 20012    72

smallnuclist <- c("U2|SNORA|SNORD|RPL|RNU|RNA5S|RPS|RNVU") 
SET1_exp_counts <- SET1_exp_counts[!grepl(smallnuclist, rownames(SET1_exp_counts)),]
dim(SET1_exp_counts) # 18923 x 72

SET2_pheno<-pheno_data[pheno_data$SET=="2",]
dim(SET2_pheno)
SET2_exp_counts<-cpm_genes[,pheno_data$SET=="2"]
dim(SET2_exp_counts)# 20012    76

smallnuclist <- c("U2|SNORA|SNORD|RPL|RNU|RNA5S|RPS|RNVU") 
SET2_exp_counts <- SET2_exp_counts[!grepl(smallnuclist, rownames(SET2_exp_counts)),]
dim(SET2_exp_counts) # 18923 x 76


## Split data into SET1 and SET2
SET1_pheno<-pheno_data[pheno_data$SET=="1",]

SET1_pheno$DX[SET1_pheno$DX=="Imp" | SET1_pheno$DX== "MS" | SET1_pheno$DX=="IMP" | SET1_pheno$DX=="Other"] = "Other"
pheno_data1 = data.frame(DIAGNOSIS=SET1_pheno$DX, SEX=SET1_pheno$SEX,RACE=SET1_pheno$RACE,AGE=SET1_pheno$AGE)

pheno_data1[1:10,1:4]


# Create AGE_GROUP before plotting
pheno_data1$AGE_GROUP <- cut(pheno_data1$AGE, 
                             breaks = c(40, 60, 70, 80, 100), 
                             labels = c("40–60", "61–70", "71–80", "81–100"), 
                             include.lowest = TRUE)

#### Supp Fig 2
# Define colors
diag_colors <- brewer.pal(length(unique(pheno_data1$DIAGNOSIS)), "Set1")
sex_colors <- c("pink", "lightblue")
race_colors <- brewer.pal(length(unique(pheno_data1$RACE)), "Pastel1")
age_colors <- brewer.pal(length(unique(pheno_data1$AGE_GROUP)), "Accent")

# Function to draw pie chart with legend on the right
create_pie_chart <- function(data, main_title, colors) {
  counts <- table(data)
  props <- round(100 * counts / sum(counts), 1)
  labels <- paste(names(counts), "(", props, "%)", sep = "")
  slice_labels <- names(counts)  # These are the labels shown on the pie
  
  par(xpd = FALSE)  # Clip to plot region
  pie(counts, main = main_title, col = colors, radius = 0.8,cex=1.3,cex.main = 1.5)
  
  par(xpd = TRUE)  # Allow legend outside
  legend("topright", inset = c(-0.05, 0), legend = labels, fill = colors, cex = 1.2, bty = "n")
}

# Open TIFF device Supp FIG 2
tiff("Plots/SET1_combined_pie_charts_legends_outside.tiff", width = 10, height = 10, units = "in", res = 300)

# 2x2 layout
par(mfrow = c(2, 2), mar = c(2, 4, 2, 5))  # Extra right margin for legend

# Draw pie charts with adjusted legend
create_pie_chart(pheno_data1$DIAGNOSIS, "Diagnosis", diag_colors)
create_pie_chart(pheno_data1$SEX, "Sex", sex_colors)
create_pie_chart(pheno_data1$RACE, "Race", race_colors)
create_pie_chart(pheno_data1$AGE_GROUP, "Age Group", age_colors)

# Close device
dev.off()


###set2 pie chart
SET2_pheno<-pheno_data[pheno_data$SET=="2",]

#SET1_pheno$DX[SET1_pheno$DX=="Imp" | SET1_pheno$DX== "MS" | SET1_pheno$DX=="IMP" | SET1_pheno$DX=="Other"] = "Other"
pheno_data2 = data.frame(DIAGNOSIS=SET2_pheno$DX, SEX=SET2_pheno$SEX,RACE=SET2_pheno$RACE,AGE=SET2_pheno$AGE)

pheno_data2[1:10,1:4]


# Create AGE_GROUP before plotting
pheno_data2$AGE_GROUP <- cut(pheno_data2$AGE, 
                             breaks = c(40, 60, 70, 80, 100), 
                             labels = c("40–60", "61–70", "71–80", "81–100"), 
                             include.lowest = TRUE)

# Load required package
# Define colors
diag_colors <- c("#8DA0CB", "#E78AC3")
sex_colors <- c("pink", "lightblue")
race_colors <- brewer.pal(length(unique(pheno_data2$RACE)), "Pastel1")
age_colors <- brewer.pal(length(unique(pheno_data2$AGE_GROUP)), "Accent")

# Function to draw pie chart with legend on the right
create_pie_chart <- function(data, main_title, colors) {
  counts <- table(data)
  props <- round(100 * counts / sum(counts), 1)
  labels <- paste(names(counts), "(", props, "%)", sep = "")
  slice_labels <- names(counts)  # These are the labels shown on the pie
  
  par(xpd = FALSE)  # Clip to plot region
  pie(counts, main = main_title, col = colors, radius = 0.8,cex=1.3,cex.main = 1.5)
  
  par(xpd = TRUE)  # Allow legend outside
  legend("topright", inset = c(-0.05, 0), legend = labels, fill = colors, cex = 1.2, bty = "n")
}

# Open TIFF device
tiff("PLots/SET2_combined_pie_charts_legends_outside.tiff", width = 10, height = 10, units = "in", res = 300)

# 2x2 layout
par(mfrow = c(2, 2), mar = c(2, 4, 2, 5))  # Extra right margin for legend

# Draw pie charts with adjusted legend
create_pie_chart(pheno_data2$DIAGNOSIS, "Diagnosis", diag_colors)
create_pie_chart(pheno_data2$SEX, "Sex", sex_colors)
create_pie_chart(pheno_data2$RACE, "Race", race_colors)
create_pie_chart(pheno_data2$AGE_GROUP, "Age Group", age_colors)
# Close device
dev.off()

################################################################################
####  STEP 2: Split the data into train and remaining sets based on the target variable (CAT)
################################################################################

set.seed(1230)
index <- createDataPartition(SET1_pheno$CAT, p = 0.8, list = FALSE)

# Split both feature and target datasets for training data
trainData_counts <- SET1_exp_counts[, index]
trainingData_pheno <- SET1_pheno$CAT[index]
pheno_train<- SET1_pheno[index,]
dim(pheno_train)

# Split the remaining data for validation and testing
remainingData_counts <- SET1_exp_counts[, -index]
remainingData_pheno <- SET1_pheno$CAT[-index]
pheno_valid <- SET1_pheno[-index,]
dim(pheno_valid)

# Split the remaining data into validation (60%) and test (40%) sets
# remaining_index <- createDataPartition(remainingData_pheno, p = 0.6, list = FALSE)

validationData_counts <- remainingData_counts#[, remaining_index]
validationData_pheno <- remainingData_pheno#[remaining_index]

# Print dimensions of splits
cat("Training data dimensions:", dim(trainData_counts), "\n")
cat("Validation data dimensions:", dim(validationData_counts), "\n")

# Names are wrong
cpm_train<-trainData_counts
cpm_valid <-validationData_counts
cpm_test <- SET2_exp_counts


################################################################################
## New figure 2
plot_pTauT<-ggplot(pheno_train, aes(x = AGE, y = (pTAU.ABETA),ymin =0.003, ymax=0.3)) +
  geom_point(aes(color = SEX),size=2.5) +
  facet_wrap(~ SEX) +
  labs(title = "pTau/Abeta42 ratios vs. Age by Sex",
       x = "Age",
       y = "pTau/Abeta42") +
  scale_y_log10(breaks=c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1), labels=c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1)) +
  geom_hline(yintercept=(0.025), linetype="solid", color="black", linewidth=1) +
  theme_minimal()+theme(plot.title = element_text(face = "bold",size=14),
                        axis.title = element_text(face = "bold", size=14),
                        axis.text = element_text(size=14),
                        strip.text = element_text(face = "bold",size=14),
                        legend.title = element_text(face = "bold",size=14),
                        legend.text = element_text(size=14))

tiff(file = "Plots/01.PTAUabeta_SET1_Train.tiff", width = 5, height = 4, units = "in", res = 300)
print(plot_pTauT)
dev.off()

plot_pTauV<-ggplot(pheno_valid, aes(x = AGE, y = (pTAU.ABETA),ymin =0.003, ymax=0.3)) +
  geom_point(aes(color = SEX),size=2.5) +
  facet_wrap(~ SEX) +
  labs(title = "pTau/Abeta42 ratios vs. Age by Sex",
       x = "Age",
       y = "pTau/Abeta42") +
  scale_y_log10(breaks=c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1), labels=c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1)) +
  geom_hline(yintercept=(0.025), linetype="solid", color="black", linewidth=1) +
  theme_minimal()+theme(plot.title = element_text(face = "bold",size=14),
                        axis.title = element_text(face = "bold", size=14),
                        axis.text = element_text(size=14),
                        strip.text = element_text(face = "bold",size=14),
                        legend.title = element_text(face = "bold",size=14),
                        legend.text = element_text(size=14))

tiff(file = "Plots/01.PTAUabeta_SET1_Valid.tiff", width = 5, height = 4, units = "in", res = 300)
print(plot_pTauV)
dev.off()

plot_pTau<-ggplot(SET2_pheno, aes(x = AGE, y = (pTAU.ABETA),ymin =0.003, ymax=0.3)) +
  geom_point(aes(color = SEX),size=2.5) +
  facet_wrap(~ SEX) +
  labs(title = "pTau/Abeta42 ratios vs. Age by Sex",
       x = "Age",
       y = "pTau/Abeta42") +
  scale_y_log10(breaks=c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1), labels=c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1)) +
  geom_hline(yintercept=(0.025), linetype="solid", color="black", linewidth=1) +
  theme_minimal()+theme(plot.title = element_text(face = "bold",size=14),
                        axis.title = element_text(face = "bold", size=14),
                        axis.text = element_text(size=14),
                        strip.text = element_text(face = "bold",size=14),
                        legend.title = element_text(face = "bold",size=14),
                        legend.text = element_text(size=14))

tiff(file = "Plots/01.PTAUabeta_SET2_Test.tiff", width = 5, height = 4, units = "in", res = 300)
print(plot_pTau)
dev.off()


################################################################################
####. STEP 3:  CREATE SIGNIFICANT GENES DATASHEET         ######################
################################################################################


# Assuming 'cpm_train' is a matrix (genes x samples) and 
# 'SET1_pheno$ptau.abeta[index]' is a vector of outcomes
correlation_results <- apply(trainData_counts, 1, function(x) {
  cor.test(x, SET1_pheno$ptau[index], method = "pearson")
})

# Extract the correlation coefficients and p-values
correlations <- sapply(correlation_results, function(res) res$estimate)
p_values <- sapply(correlation_results, function(res) res$p.value)

adjusted_p_values_bh <- p.adjust(p_values, method = "BH")

results <- data.frame(
  correlation = correlations,
  p_value = p_values,
  p_value_adj_bh = adjusted_p_values_bh
)

# Check the results
head(results)
results
summary(results$correlation)

# subset to absolute correlation > 0.6
results[1:4,]
corgenes<-results[abs(results$correlation)>0.50 ,]
corgenes2<-corgenes[complete.cases(corgenes),]
corgenes2$correlation
corgenes2$genes <- sub("\\.cor$", "", rownames(corgenes2))
corgenes2 #check genes
dim(corgenes2)

correlation<-data.frame(cor=corgenes2$correlation,genes=corgenes2$genes)
genes<-correlation$genes
write.csv(corgenes2,"onlyPTAU_SET1traingenes.csv")
write.csv(SET1_pheno,"SET1_pheno.csv")


################################################################################

## Subset train data based on significant genes
## if imported list use this
#genes<-sig_genes$genes

train <- cpm_train[rownames(cpm_train) %in% genes,]
dim(train)
# Subset validation data based on filtered genes
validation <- cpm_valid[rownames(cpm_train) %in% genes, ]
dim(validation)
str(validation)

testSet2<-cpm_test[rownames(cpm_train) %in% genes, ]
str(testSet2)

################################################################################
#### STEP 4: Center and scale the training data (do this only on the training set)
################################################################################
set.seed(1230)
preProcValues <- preProcess(t(train), method= c("center", "scale"))

# Apply the centering and scaling transformation to training, validation, and test sets
train1 <- predict(preProcValues, t(train))
validation1 <- predict(preProcValues, t(validation))
#test1 <- predict(preProcValues, t(test))
test_SET2<-predict(preProcValues,t(testSet2))

# Check the dimensions to make sure they are correct
cat("Training set dimensions after centering and scaling:", dim(train1), "\n")
cat("Validation set dimensions after centering and scaling:", dim(validation1), "\n")
cat("Test set dimensions after centering and scaling:", dim(test_SET2), "\n")

# Prepare the final training, validation, and test data by adding the target variable (CAT)
# Add target variable (CAT) for training data
group_train <- trainingData_pheno  # Training group (CAT)
train <- data.frame(train1, CATs = group_train)  # Add CATs to the training data

# Add target variable (CAT) for validation data
group_validation <- validationData_pheno # Validation group (CAT)
valid <- data.frame(validation1, CATs = group_validation)  # Add CATs to the validation data

# Add target variable (CAT) for test data
testSET2<-data.frame(test_SET2,CATs=SET2_pheno$CAT)

# Check final dimensions of the training, validation, and test sets
cat("Final training set dimensions:", dim(train), "\n")
cat("Final validation set dimensions:", dim(validation), "\n")
cat("SET2 dimensions:", dim(testSET2), "\n")

################################################################################
##### STEP 5: TRAINING MODELS ON THE TRAIN DATA        #########################
################################################################################

# Create a train control object
fitControl <- trainControl(
  method = 'cv',                   # k-fold cross validation
  number = 10,                   # number of folds
  savePredictions = 'final',    # saves predictions for optimal tuning parameter
  classProbs = TRUE, verboseIter = FALSE, # should class probabilities returned
  summaryFunction = twoClassSummary)  # results summary function

dim(train); dim(valid)
# Now run various models using the tuning grid

################################################################################
## A: Random Forest (RF)
set.seed(1234)  
fit_rf <- train(CATs ~ ., data = train, method = "rf", tuneLength=20, trControl = fitControl,metric="ROC")
print(fit_rf)
summary(fit_rf)

## B: Ranger
set.seed(1234)
fit_ranger <- train(CATs ~ ., data = train, importance="impurity", method = "ranger", tuneLength=10, trControl = fitControl,metric="ROC")
# Print the results
print(fit_ranger)
summary(fit_ranger)

## C: Neural-NET
set.seed(1234)
# Example with specific tuning for the neural network
fit_nnet <- train(CATs ~ ., data = train,
                  method = "nnet", 
                  trControl = fitControl, 
                  tuneGrid = expand.grid(.size = c(5, 8, 10, 12, 15), 
                                         .decay = c(0.001, 0.01, 0.1)),
                  linout = FALSE,
                  trace = FALSE,metric="ROC")
print(fit_nnet)

## D: GLM-NET
grids = expand.grid(alpha = seq(0.1, 1, 0.01), lambda = seq(0.01, 0.1, length = 100)) 
set.seed(1234)
fit_glm <- train(CATs ~ ., data = train, method = "glmnet", tuneGrids=grids, trControl = fitControl,metric="ROC")
print(fit_glm)

## E: Partial Least Squares (PLS)
set.seed(1234)
fit_pls <- train(CATs ~., data = train, method = "pls", tuneLength=10, trControl = fitControl,metric="ROC")
print(fit_pls)

## F: Support Vector machine (Radial) SVMRADIAL
set.seed(1234)
fit_svm <- train(CATs ~., data = train, method = "svmRadial", tuneLength=20, trControl = fitControl,metric="ROC")
print (fit_svm)

## G: RPART
set.seed(1234)
fit_rpart <- train(CATs ~., data = train, method = "rpart", tuneLength=10, trControl = fitControl,metric="ROC")
print(fit_rpart)

## H: MARS
set.seed(1234)
fit_earth <- train(CATs ~., data = train, method = "earth", tuneLength=10, trControl = fitControl,metric="ROC")
print(fit_earth)


## I: Gradient Boost
set.seed(1234) 
fit_GBM = train(CATs~ ., data=train, method='gbm', tuneLength=10, trControl = fitControl)
summary(fit_GBM)
print(fit_GBM)
################################################################################


################################################################################
#######      STEP.6. Compare Models             ################################
################################################################################

## Assess training ROC using MLeval Fig4 ROC training SET1
tiff("Plots/roc_ModelsSET1.tiff",width=4,height=4, unit= "in", res = 300 )
res <- evalm(list(fit_rf,fit_glm,fit_earth,
                  fit_ranger,fit_pls,fit_nnet,fit_rpart,fit_svm, fit_GBM),
             gnames=c("rf","glmnet","earth",
                      "ranger","pls","nnet","rpart","svm","GBM"))
dev.off()

##  Then compare model performances using resample()
models_compare <- resamples(list(rf=fit_rf,glm=fit_glm,earth=fit_earth,
                                 ranger=fit_ranger,pls=fit_pls,nnet=fit_nnet,rpart=fit_rpart,svm=fit_svm, GBM=fit_GBM))

# Summary of the models performances
sink("RESULTS/Models_Compare_Set1_Train.txt")
print(summary(models_compare))
sink()

#tiff(file="SET1_train_Box_comp.tiff", unit= "in", res = 500, width = 9, height = 6)
bwplot(models_compare)
#dev.off()

# Modify font size and make text bold
par.settings <- list(
  axis.text = list(cex = 1.5, font = 1),  # Increase axis text size and make it bold
  axis.title = list(cex = 1.5, font = 2),  # Increase axis title size and make it bold
  strip.text = list(cex = 1.5, font = 2),  # Increase strip text size and make it bold
  layout.heights = list(top.padding = 10),  # Adjust the top padding for space if needed
  layout.widths = list(left.padding = 10)  # Adjust left padding if needed
)
scales <- list(x=list(relation="free"), y=list(relation="free"))

# Plot with the modified settings
tiff(file="Plots/SET1_train_Box_comp.tiff", unit="in", res=300, width=9, height=6)
bwplot(models_compare, scales=scales, par.settings=par.settings)
dev.off()

################################################################################
### STEP.7 Testing the models
################################################################################
## A. Test RF model
# test on train data-categorical
sink("RESULTS/RF_train_valid.txt")
output_train_rf<- predict(fit_rf ,train)
cm_train_rf <- confusionMatrix(reference = as.factor(train$CATs), data =  as.factor(output_train_rf),mode='everything', positive='YES')
# Print the results
print(cm_train_rf)

##predict on validation(0.2)
output_valid_rf  <- predict(fit_rf, valid)
cm_valid_rf <- confusionMatrix(reference = as.factor(valid$CATs), data = as.factor(output_valid_rf), mode='everything', positive='YES') 
print(cm_valid_rf)

sink()

# calculate probabilities and ROC using Probs
output_train_SET1_rf_P  <- predict(fit_rf, train, type="prob")
roc_rf <- roc(train$CATs, output_train_SET1_rf_P$YES)
plot(roc_rf)
roc_rf


output_valid_SET1_rf_P  <- predict(fit_rf, valid, type="prob")
roc_rfV <- roc(valid$CATs, output_valid_SET1_rf_P$YES)
plot(roc_rfV)
roc_rfV


## B. Test Ranger model
# test on train data-categorical
#predict on train(0.5)
output_train_ranger<- predict(fit_ranger ,train)
cm_train_ranger <- confusionMatrix(reference = as.factor(train$CATs), data =  as.factor(output_train_ranger),mode='everything', positive='YES')
print(cm_train_ranger)

# calculate probabilities and train ROC using Probs
output_train_SET1_ranger_P  <- predict(fit_ranger, train, type="prob")
str(output_train_SET1_ranger_P)
roc_ranger <- roc(train$CATs, output_train_SET1_ranger_P$YES)
plot(roc_ranger)
roc_ranger

##predict on validation(0.3)
output_valid_ranger  <- predict(fit_ranger, valid)
cm_valid_ranger <- confusionMatrix(reference = as.factor(valid$CATs), data = as.factor(output_valid_ranger), mode='everything', positive='YES') 
print(cm_valid_ranger)


## C. Test NNET model
# test on train data-categorical
output_train_nnet<- predict(fit_nnet ,train)
cm_train_nnet <- confusionMatrix(reference = as.factor(train$CATs), data =  as.factor(output_train_nnet),mode='everything', positive='YES')
print(cm_train_nnet)

# calculate probabilities and train ROC using Probs
output_train_SET1_nnet_P  <- predict(fit_nnet, train, type="prob")
str(output_train_SET1_nnet_P)
roc_nnet <- roc(train$CATs, output_train_SET1_nnet_P$YES)
plot(roc_nnet)
roc_nnet

##predict on validation(0.3)
output_valid_nnet  <- predict(fit_nnet, valid)
cm_valid_nnet <- confusionMatrix(reference = as.factor(valid$CATs), data = as.factor(output_valid_nnet), mode='everything', positive='YES') 
print(cm_valid_nnet)


## D. Test GLM model
# test on train data-categorical
output_train_glm<- predict(fit_glm ,train)
cm_train_glm <- confusionMatrix(reference = as.factor(train$CATs), data =  as.factor(output_train_glm),mode='everything', positive='YES')
print(cm_train_glm)

# calculate probabilities and train ROC using Probs
output_train_SET1_glm_P  <- predict(fit_glm, train, type="prob")
str(output_train_SET1_glm_P)
roc_glm <- roc(train$CATs, output_train_SET1_glm_P$YES)
plot(roc_glm)
roc_glm

##predict on validation(0.2)
output_valid_glm  <- predict(fit_glm, valid)
cm_valid_glm <- confusionMatrix(reference = as.factor(valid$CATs), data = as.factor(output_valid_glm), mode='everything', positive='YES') 
print(cm_valid_glm)


## E. Test PLS model
# test on train data-categorical
output_train_pls<- predict(fit_pls ,train)
cm_train_pls <- confusionMatrix(reference = as.factor(train$CATs), data =  as.factor(output_train_pls),mode='everything', positive='YES')
print(cm_train_pls)

# calculate probabilities and train ROC using Probs
output_train_SET1_pls_P  <- predict(fit_pls, train, type="prob")
str(output_train_SET1_pls_P)
roc_pls <- roc(train$CATs, output_train_SET1_pls_P$YES)
plot(roc_pls)
roc_pls

##predict on validation(0.2)
output_valid_pls  <- predict(fit_pls, valid)
cm_valid_pls <- confusionMatrix(reference = as.factor(valid$CATs), data = as.factor(output_valid_pls), mode='everything', positive='YES') 
print(cm_valid_pls)


## F. Test SVM-Radial model
# test on train data-categorical
output_train_svm<- predict(fit_svm ,train)
cm_train_svm <- confusionMatrix(reference = as.factor(train$CATs), data =  as.factor(output_train_svm),mode='everything', positive='YES')
print(cm_train_svm)

# calculate probabilities and train ROC using Probs
output_train_SET1_svm_P  <- predict(fit_svm, train, type="prob")
str(output_train_SET1_svm_P)
roc_svm <- roc(train$CATs, output_train_SET1_svm_P$YES)
plot(roc_svm)
roc_svm

##predict on validation(0.2)
output_valid_svm  <- predict(fit_svm, valid)
cm_valid_svm <- confusionMatrix(reference = as.factor(valid$CATs), data = as.factor(output_valid_svm), mode='everything', positive='YES') 
print(cm_valid_svm)


## G. Test RPART model
# test on train data-categorical
output_train_rpart<- predict(fit_rpart ,train)
cm_train_rpart <- confusionMatrix(reference = as.factor(train$CATs), data =  as.factor(output_train_rpart),mode='everything', positive='YES')
print(cm_train_rpart)

# calculate probabilities and train ROC using Probs
output_train_SET1_rpart_P  <- predict(fit_rpart, train, type="prob")
str(output_train_SET1_rpart_P)
roc_rpart <- roc(train$CATs, output_train_SET1_rpart_P$YES)
plot(roc_rpart)
roc_rpart

##predict on validation(0.2)
output_valid_rpart  <- predict(fit_rpart, valid)
cm_valid_rpart <- confusionMatrix(reference = as.factor(valid$CATs), data = as.factor(output_valid_rpart), mode='everything', positive='YES') 
print(cm_valid_rpart)


## H. Test Mars model
# test on train data-categorical
# test on train 0.5
output_train_earth<- predict(fit_earth ,train)
cm_train_earth <- confusionMatrix(reference = as.factor(train$CATs), data =  as.factor(output_train_earth),mode='everything', positive='YES')
print(cm_train_earth)

# calculate probabilities and create ROC
output_train_SET1_earth_P  <- predict(fit_earth, train, type="prob")
str(output_train_SET1_earth_P)
roc_earth <- roc(train$CATs, output_train_SET1_earth_P$YES)
plot(roc_earth)
roc_earth

##predict on validation(0.2)
output_valid_earth  <- predict(fit_earth, valid)
cm_valid_earth <- confusionMatrix(reference = as.factor(valid$CATs), data = as.factor(output_valid_earth), mode='everything', positive='YES') 
print(cm_valid_earth)


## I. # Run GBM model to predict train 0.5
sink("RESULTS/GBM_train_valid.txt")
output_train_GBM<- predict(fit_GBM ,train)
cm_train_GBM <- confusionMatrix(reference = as.factor(train$CATs), data =  as.factor(output_train_GBM),mode='everything', positive='YES')
print(cm_train_GBM)

##predict on validation(0.2)
output_valid_GBM  <- predict(fit_GBM, valid)
cm_valid_GBM <- confusionMatrix(reference = as.factor(valid$CATs), data = as.factor(output_valid_GBM), mode='everything', positive='YES') 
print(cm_valid_GBM)
sink()

# calculate probabilities and create ROC
output_train_SET1_GBM_P  <- predict(fit_GBM, train, type="prob")
str(output_train_SET1_GBM_P)
roc_GBM <- roc(train$CATs, output_train_SET1_GBM_P$YES)
plot(roc_GBM)
roc_GBM

output_valid_SET1_GBM_P  <- predict(fit_GBM, valid, type="prob")
roc_GBMV <- roc(valid$CATs, output_valid_SET1_GBM_P$YES)
plot(roc_GBMV)
roc_GBMV


################################################################################
###########  STEP 8 Testing on SET 2 DATA             ##########################
################################################################################
## The RF model was best in training, but we add all models for completeness. 

##now test rf Model on SET2
output_test_SET2_rf  <- predict(fit_rf, testSET2)
cm_test_SET2_rf <- confusionMatrix(reference = as.factor(testSET2$CATs), data = as.factor(output_test_SET2_rf), mode='everything', positive='YES') 

#get pvalue for ROC-AUC
output_test_SET2_rf$OBS <-  ifelse(testSET2$CAT == "YES", 1, 0)
roc_pvalue_rf <- roc.area(output_test_SET2_rf$OBS,as.numeric(output_test_SET2_rf$YES))

sink("RESULTS/RF_Test.txt")
print(cm_test_SET2_rf)
cat("ROC-AUC analysis\n")
print(roc_pvalue_rf)
sink()


##now test ranger Model on SET2
output_test_SET2_ranger  <- predict(fit_ranger, testSET2)
cm_test_SET2_ranger <- confusionMatrix(reference = as.factor(testSET2$CATs), data = as.factor(output_test_SET2_ranger), mode='everything', positive='YES') 
print(cm_test_SET2_ranger)

##now test nnet Model on SET2
output_test_SET2_nnet  <- predict(fit_nnet, testSET2)
cm_test_SET2_nnet <- confusionMatrix(reference = as.factor(testSET2$CATs), data = as.factor(output_test_SET2_nnet), mode='everything', positive='YES') 
print(cm_test_SET2_nnet)

##now test pls Model on SET2
output_test_SET2_pls  <- predict(fit_pls, testSET2)
cm_test_SET2_pls <- confusionMatrix(reference = as.factor(testSET2$CATs), data = as.factor(output_test_SET2_pls), mode='everything', positive='YES') 
print(cm_test_SET2_pls)

##now test svm Model on SET2
output_test_SET2_svm  <- predict(fit_svm, testSET2)
cm_test_SET2_svm <- confusionMatrix(reference = as.factor(testSET2$CATs), data = as.factor(output_test_SET2_svm), mode='everything', positive='YES') 
print(cm_test_SET2_svm)

##now test rpart Model on SET2
output_test_SET2_rpart  <- predict(fit_rpart, testSET2)
cm_test_SET2_rpart <- confusionMatrix(reference = as.factor(testSET2$CATs), data = as.factor(output_test_SET2_rpart), mode='everything', positive='YES') 
print(cm_test_SET2_rpart)

##now test earth Model on SET2
output_test_SET2_earth  <- predict(fit_earth, testSET2)
cm_test_SET2_earth <- confusionMatrix(reference = as.factor(testSET2$CATs), data = as.factor(output_test_SET2_earth), mode='everything', positive='YES') 
print(cm_test_SET2_earth)

##now test GBM Model on SET2
output_test_SET2_GBM  <- predict(fit_GBM, testSET2)
cm_test_SET2_GBM <- confusionMatrix(reference = as.factor(testSET2$CATs), data = as.factor(output_test_SET2_GBM), mode='everything', positive='YES') 

sink("RESULTS/GBM_Test.txt")
print(cm_test_SET2_GBM)
sink()


################################################################################
###############.    STEP 9 Tuning the probabilities to improve performance  ####
################################################################################
## AUC used a Prob of 0.5 to call samples, but is this the best cutoff to use?  
# we can re-train the cutoff using AUC :)

### A: Random Forest
best_threshold <- coords(roc_rf, "best", ret="threshold")
print(best_threshold)

#best=0.529
## USE this to re-code probabilities

# calculate probabilities using Probs apply new threshold
output_valid_rfP<- predict(fit_rf ,valid, type="prob")
output_valid_rfP$pCAT <- ifelse(output_valid_rfP$YES > 0.53, "YES", "NO") 
output_valid_rfP<- data.frame(output_valid_rfP)
cm_valid_rfP <- confusionMatrix(reference = as.factor(valid$CATs), data =  as.factor(output_valid_rfP[,3]),mode='everything', positive='YES')
print(cm_valid_rfP)
#draw ROC
roc_rfv <- roc(valid$CATs, output_valid_rfP$YES)
plot(roc_rfv)
roc_rfv

# calculate probabilities using Probs apply new threshold
output_test_SET2_rf_P  <- predict(fit_rf, testSET2, type="prob")
output_test_SET2_rf_P$pCAT <- ifelse(output_test_SET2_rf_P$YES > 0.53, "YES", "NO") 
output_test_SET2_rf_P<- data.frame(output_test_SET2_rf_P)
cm_test_SET2_rf_P <- confusionMatrix(reference = as.factor(testSET2$CATs), data = as.factor(output_test_SET2_rf_P[,3]), mode='everything', positive='YES') 
print(cm_test_SET2_rf_P)

#get pvalue for ROC-AUC
output_test_SET2_rf_P$OBS <-  ifelse(testSET2$CAT== "YES", 1, 0)
roc_pvalue_rf_P <- roc.area(output_test_SET2_rf_P$OBS,as.numeric(output_test_SET2_rf_P$YES))

sink("RESULTS/RF_tuned_Test.txt")
print(cm_test_SET2_rf_P)
cat("ROC-AUC analysis\n")
print(roc_pvalue_rf_P)
sink()

#draw ROC
roc_rfP <- roc(testSET2$CATs, output_test_SET2_rf_P$YES, ci=TRUE)
plot(roc_rfP)
roc_rfP

tiff(file="PLOTS/Tuned_RF_Test_ROC.tiff", unit="in", res=500, width=3, height=3)
plot(roc_rfP)
dev.off()

### B: Ranger
best_threshold <- coords(roc_ranger, "best", ret="threshold")
print(best_threshold)
#best=0.56

# calculate probabilities and train ROC using Probs
output_valid_SET1_ranger_P  <- predict(fit_ranger, valid, type="prob")
output_valid_SET1_ranger_P$pCAT <- ifelse(output_valid_SET1_ranger_P$YES > 0.56, "YES", "NO") 
output_valid_SET1_ranger_P<- data.frame(output_valid_SET1_ranger_P)
cm_valid_SET1_ranger_P <- confusionMatrix(reference = as.factor(valid$CATs), data =  as.factor(output_valid_SET1_ranger_P[,3]),mode='everything', positive='YES')
print(cm_valid_SET1_ranger_P)

#draw ROC
roc_rangerv <- roc(valid$CATs, output_valid_SET1_ranger_P$YES)
plot(roc_rangerv)
roc_rangerv

# calculate probabilities using Probs apply new threshold
output_test_SET2_ranger_P  <- predict(fit_ranger, testSET2, type="prob")
output_test_SET2_ranger_P$pCAT <- ifelse(output_test_SET2_ranger_P$YES > 0.56, "YES", "NO") 
output_test_SET2_ranger_P<- data.frame(output_test_SET2_ranger_P)
cm_test_SET2_ranger_P <- confusionMatrix(reference = as.factor(testSET2$CATs), data = as.factor(output_test_SET2_ranger_P[,3]), mode='everything', positive='YES') 
print(cm_test_SET2_ranger_P)

#draw ROC
roc_rangerP <- roc(testSET2$CATs, output_test_SET2_ranger_P$YES)
plot(roc_rangerP)
roc_rangerP


### C: NNET
best_threshold <- coords(roc_nnet, "best", ret="threshold")
print(best_threshold)
#best=0.4996

# calculate probabilities and train ROC using Probs
output_valid_SET1_nnet_P  <- predict(fit_nnet, valid, type="prob")
output_valid_SET1_nnet_P$pCAT <- ifelse(output_valid_SET1_nnet_P$YES > 0.4996, "YES", "NO") 
output_valid_SET1_nnet_P<- data.frame(output_valid_SET1_nnet_P)
cm_valid_SET1_nnet_P <- confusionMatrix(reference = as.factor(valid$CATs), data =  as.factor(output_valid_SET1_nnet_P[,3]),mode='everything', positive='YES')
print(cm_valid_SET1_nnet_P)

#draw ROC
roc_nnetv <- roc(valid$CATs, output_valid_SET1_nnet_P$YES)
plot(roc_nnetv)
roc_nnetv

# calculate probabilities using Probs apply new threshold
output_test_SET2_nnet_P  <- predict(fit_nnet, testSET2, type="prob")
output_test_SET2_nnet_P$pCAT <- ifelse(output_test_SET2_nnet_P$YES > 0.4996, "YES", "NO") 
output_test_SET2_nnet_P<- data.frame(output_test_SET2_nnet_P)
cm_test_SET2_nnet_P <- confusionMatrix(reference = as.factor(testSET2$CATs), data = as.factor(output_test_SET2_nnet_P[,3]), mode='everything', positive='YES') 
print(cm_test_SET2_nnet_P)

#draw ROC
roc_nnetP <- roc(testSET2$CATs, output_test_SET2_nnet_P$YES)
plot(roc_nnetP)
roc_nnetP

### D: GLM
best_threshold <- coords(roc_glm, "best", ret="threshold")
print(best_threshold)
#best=0.6549

# calculate probabilities and train ROC using Probs
output_valid_SET1_glm_P  <- predict(fit_glm, valid, type="prob")
output_valid_SET1_glm_P$pCAT <- ifelse(output_valid_SET1_glm_P$YES > 0.6549, "YES", "NO") 
output_valid_SET1_glm_P<- data.frame(output_valid_SET1_glm_P)
cm_valid_SET1_glm_P <- confusionMatrix(reference = as.factor(valid$CATs), data =  as.factor(output_valid_SET1_glm_P[,3]),mode='everything', positive='YES')
print(cm_valid_SET1_glm_P)

#draw ROC
roc_glmv <- roc(valid$CATs, output_valid_SET1_glm_P$YES)
plot(roc_glmv)
roc_glmv

##now test glm Model on SET2
output_test_SET2_glm  <- predict(fit_glm, testSET2)
cm_test_SET2_glm <- confusionMatrix(reference = as.factor(testSET2$CATs), data = as.factor(output_test_SET2_glm), mode='everything', positive='YES') 
print(cm_test_SET2_glm)

# calculate probabilities using Probs apply new threshold
output_test_SET2_glm_P  <- predict(fit_glm, testSET2, type="prob")
output_test_SET2_glm_P$pCAT <- ifelse(output_test_SET2_glm_P$YES > 0.6549, "YES", "NO") 
output_test_SET2_glm_P<- data.frame(output_test_SET2_glm_P)
cm_test_SET2_glm_P <- confusionMatrix(reference = as.factor(testSET2$CATs), data = as.factor(output_test_SET2_glm_P[,3]), mode='everything', positive='YES') 
print(cm_test_SET2_glm_P)

#draw ROC
roc_glmP <- roc(testSET2$CATs, output_test_SET2_glm_P$YES)
plot(roc_glmP)
roc_glmP


### E: PLS
best_threshold <- coords(roc_pls, "best", ret="threshold")
print(best_threshold)
#best=0.5408

# calculate probabilities and train ROC using Probs
output_valid_SET1_pls_P  <- predict(fit_pls, valid, type="prob")
output_valid_SET1_pls_P$pCAT <- ifelse(output_valid_SET1_pls_P$YES > 0.5408, "YES", "NO") 
output_valid_SET1_pls_P<- data.frame(output_valid_SET1_pls_P)
cm_valid_SET1_pls_P <- confusionMatrix(reference = as.factor(valid$CATs), data =  as.factor(output_valid_SET1_pls_P[,3]),mode='everything', positive='YES')
print(cm_valid_SET1_pls_P)

#draw ROC
roc_plsv <- roc(valid$CATs, output_valid_SET1_pls_P$YES)
plot(roc_plsv)
roc_plsv

# calculate probabilities using Probs apply new threshold
output_test_SET2_pls_P  <- predict(fit_pls, testSET2, type="prob")
output_test_SET2_pls_P$pCAT <- ifelse(output_test_SET2_pls_P$YES > 0.5408, "YES", "NO") 
output_test_SET2_pls_P<- data.frame(output_test_SET2_pls_P)
cm_test_SET2_pls_P <- confusionMatrix(reference = as.factor(testSET2$CATs), data = as.factor(output_test_SET2_pls_P[,3]), mode='everything', positive='YES') 
print(cm_test_SET2_pls_P)

#draw ROC
roc_plsP <- roc(testSET2$CATs, output_test_SET2_pls_P$YES)
plot(roc_plsP)
roc_plsP


### F: Model SVM-R
best_threshold <- coords(roc_svm, "best", ret="threshold")
print(best_threshold)
#best=0.4838

# calculate probabilities and train ROC using Probs
output_valid_SET1_svm_P  <- predict(fit_svm, valid, type="prob")
output_valid_SET1_svm_P$pCAT <- ifelse(output_valid_SET1_svm_P$YES > 0.4838, "YES", "NO") 
output_valid_SET1_svm_P<- data.frame(output_valid_SET1_svm_P)
cm_valid_SET1_svm_P <- confusionMatrix(reference = as.factor(valid$CATs), data =  as.factor(output_valid_SET1_svm_P[,3]),mode='everything', positive='YES')
print(cm_valid_SET1_svm_P)

#draw ROC
roc_svmv <- roc(valid$CATs, output_valid_SET1_svm_P$YES)
plot(roc_svmv)
roc_svmv

# calculate probabilities using Probs apply new threshold
output_test_SET2_svm_P  <- predict(fit_svm, testSET2, type="prob")
output_test_SET2_svm_P$pCAT <- ifelse(output_test_SET2_svm_P$YES > 0.4838, "YES", "NO") 
output_test_SET2_svm_P<- data.frame(output_test_SET2_svm_P)
cm_test_SET2_svm_P <- confusionMatrix(reference = as.factor(testSET2$CATs), data = as.factor(output_test_SET2_svm_P[,3]), mode='everything', positive='YES') 
print(cm_test_SET2_svm_P)

#draw ROC
roc_svmP <- roc(testSET2$CATs, output_test_SET2_svm_P$YES)
plot(roc_svmP)
roc_svmP


### G: Model RPART
# calculate probabilities and train ROC using Probs
best_threshold <- coords(roc_rpart, "best", ret="threshold")
print(best_threshold)
#best=0.6406

# calculate probabilities and train ROC using Probs
output_valid_SET1_rpart_P  <- predict(fit_rpart, valid, type="prob")
output_valid_SET1_rpart_P$pCAT <- ifelse(output_valid_SET1_rpart_P$YES > 0.6406, "YES", "NO") 
output_valid_SET1_rpart_P<- data.frame(output_valid_SET1_rpart_P)
cm_valid_SET1_rpart_P <- confusionMatrix(reference = as.factor(valid$CATs), data =  as.factor(output_valid_SET1_rpart_P[,3]),mode='everything', positive='YES')
print(cm_valid_SET1_rpart_P)

#draw ROC
roc_rpartv <- roc(valid$CATs, output_valid_SET1_rpart_P$YES)
plot(roc_rpartv)
roc_rpartv

# calculate probabilities using Probs apply new threshold
output_test_SET2_rpart_P  <- predict(fit_rpart, testSET2, type="prob")
output_test_SET2_rpart_P$pCAT <- ifelse(output_test_SET2_rpart_P$YES > 0.6406, "YES", "NO") 
output_test_SET2_rpart_P<- data.frame(output_test_SET2_rpart_P)
cm_test_SET2_rpart_P <- confusionMatrix(reference = as.factor(testSET2$CATs), data = as.factor(output_test_SET2_rpart_P[,3]), mode='everything', positive='YES') 
print(cm_test_SET2_rpart_P)

#draw ROC
roc_rpartP <- roc(testSET2$CATs, output_test_SET2_rpart_P$YES)
plot(roc_rpartP)
roc_rpartP


### H:  Model MARS
# calculate probabilities and train ROC using Probs
best_threshold <- coords(roc_earth, "best", ret="threshold")
print(best_threshold)
#best=0.5

# calculate probabilities and train ROC using Probs
output_valid_SET1_earth_P  <- predict(fit_earth, valid, type="prob")
output_valid_SET1_earth_P$pCAT <- ifelse(output_valid_SET1_earth_P$YES > 0.5, "YES", "NO") 
output_valid_SET1_earth_P<- data.frame(output_valid_SET1_earth_P)
cm_valid_SET1_earth_P <- confusionMatrix(reference = as.factor(valid$CATs), data =  as.factor(output_valid_SET1_earth_P[,3]),mode='everything', positive='YES')
print(cm_valid_SET1_earth_P)

#draw ROC
roc_earthv <- roc(valid$CATs, output_valid_SET1_earth_P$YES)
plot(roc_earthv)
roc_earthv

# calculate probabilities using Probs apply new threshold
output_test_SET2_earth_P  <- predict(fit_earth, testSET2, type="prob")
output_test_SET2_earth_P$pCAT <- ifelse(output_test_SET2_earth_P$YES > 0.5, "YES", "NO") 
output_test_SET2_earth_P<- data.frame(output_test_SET2_earth_P)
cm_test_SET2_earth_P <- confusionMatrix(reference = as.factor(testSET2$CATs), data = as.factor(output_test_SET2_earth_P[,3]), mode='everything', positive='YES') 
print(cm_test_SET2_earth_P)

#draw ROC
roc_earthP <- roc(testSET2$CATs, output_test_SET2_earth_P$YES)
plot(roc_earthP)
roc_earthP


### I:  Model GBM
# calculate probabilities and train ROC using Probs
best_threshold <- coords(roc_GBM, "best", ret="threshold")
print(best_threshold)
#best=0.499

# calculate probabilities and train ROC using Probs
output_valid_SET1_GBM_P  <- predict(fit_GBM, valid, type="prob")
output_valid_SET1_GBM_P$pCAT <- ifelse(output_valid_SET1_GBM_P$YES > 0.587, "YES", "NO") 
output_valid_SET1_GBM_P<- data.frame(output_valid_SET1_GBM_P)
cm_valid_SET1_GBM_P <- confusionMatrix(reference = as.factor(valid$CATs), data =  as.factor(output_valid_SET1_GBM_P[,3]),mode='everything', positive='YES')
print(cm_valid_SET1_GBM_P)
#draw ROC
roc_gbmv <- roc(valid$CATs, output_valid_SET1_GBM_P$YES)
plot(roc_gbmv)
roc_gbmv

# calculate probabilities using Probs apply new threshold
output_test_SET2_GBM_P  <- predict(fit_GBM, testSET2, type="prob")
output_test_SET2_GBM_P$pCAT <- ifelse(output_test_SET2_GBM_P$YES > 0.587, "YES", "NO") 
output_test_SET2_GBM_P<- data.frame(output_test_SET2_GBM_P)
cm_test_SET2_GBM_P <- confusionMatrix(reference = as.factor(testSET2$CATs), data = as.factor(output_test_SET2_GBM_P[,3]), mode='everything', positive='YES') 
print(cm_test_SET2_GBM_P)

#get pvalue for ROC-AUC
output_test_SET2_GBM_P$OBS <-  ifelse(testSET2$CAT== "YES", 1, 0)
roc_pvalue_GBM_P <- roc.area(output_test_SET2_GBM_P$OBS,as.numeric(output_test_SET2_GBM_P$YES))

sink("RESULTS/GBM_tuned_Test.txt")
print(cm_test_SET2_GBM_P)
cat("ROC-AUC analysis\n")
print(roc_pvalue_GBM_P)
sink()

#draw ROC
roc_gbmP <- roc(testSET2$CATs, output_test_SET2_GBM_P$YES, ci=TRUE)
plot(roc_gbmP)
roc_gbmP


#######FIGURES#############
#######Figure 3

######Figure 4D-E

rf_probs <- predict(fit_rf, newdata = testSET2, type = "prob")[, "YES"]
gbm_probs <- predict(fit_GBM, newdata = testSET2, type = "prob")[, "YES"]


binary_actuals <- ifelse(testSET2$CATs == "YES", 1, 0)
predicted_probs <- as.numeric(rf_probs)


# Step 1: Compute ROC with 95% CI
roc_rf <- roc(binary_actuals, predicted_probs, ci = TRUE)
auc_ci <- ci.auc(roc_rf)

# Step 2: Get optimal threshold and classification
opt_threshold <- coords(roc_rf, "best", ret = "threshold")
predicted_class <- ifelse(predicted_probs >= opt_threshold[1,], 1, 0)

# Step 3: Confusion matrix and performance metrics
conf_mat <- table(Predicted = predicted_class, Actual = binary_actuals)
TN <- conf_mat["0", "0"]
TP <- conf_mat["1", "1"]
FN <- conf_mat["0", "1"]
FP <- conf_mat["1", "0"]

accuracy <- (TP + TN) / sum(conf_mat)
sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)

# Step 4: Plot ROC curve with annotations
tiff("PLOTS/ROC_rfSET2_probCurve.tiff", width = 3, height = 3, units = "in", res = 500)  
plot(roc_rf, 
     col = "black", 
     lwd = 2, 
     legacy.axes = TRUE,
     xlab = "1 - Specificity", 
     ylab = "Sensitivity",
     main = paste0("RF Prediction ROC Curve"))

#abline(a = 0, b = 1, lty = 2, col = "gray")
# Step 5: Add text to plot
text(0.5, 0.84, labels = paste0("AUC = 0.9511"),cex=0.8)
dev.off()              

       
#########################################################################################
###########Ensembling models rf and gbm  ################################################
# Ensemble: average the probabilities
## on train
rf_probs_t <- predict(fit_rf, newdata = train, type = "prob")[, "YES"]
gbm_probs_t <- predict(fit_GBM, newdata = train, type = "prob")[, "YES"]

ensemble_probs_t <- (rf_probs_t + gbm_probs_t) / 2

roc_obj_t <- roc(response = train$CATs, predictor = ensemble_probs_t, levels = c("NO", "YES"))

# Step 2: Find optimal threshold using Youden's J (maximizes sensitivity + specificity - 1)
coords_result_t <- coords(roc_obj_t, x = "best", best.method = "youden", transpose = FALSE)

# View the threshold
coords_result_t["threshold"]
#0.514
# Final prediction
ensemble_pred_t <- ifelse(ensemble_probs_t > 0.514, "YES", "NO")


# Evaluate performance
confusionMatrix(as.factor(ensemble_pred_t), as.factor(train$CATs))

###################
## on valid
rf_probs_v <- predict(fit_rf, newdata = valid, type = "prob")[, "YES"]
gbm_probs_v <- predict(fit_GBM, newdata = valid, type = "prob")[, "YES"]

ensemble_probs_v <- (rf_probs_v + gbm_probs_v) / 2

roc_obj_v <- roc(response = valid$CATs, predictor = ensemble_probs_v, levels = c("NO", "YES"))

# Step 2: Find optimal threshold using Youden's J (maximizes sensitivity + specificity - 1)
coords_result_v <- coords(roc_obj_v, x = "best", best.method = "youden", transpose = FALSE)

# View the threshold
coords_result_v["threshold"]
#0.5917
# Final prediction
ensemble_pred_v <- ifelse(ensemble_probs_v > 0.5917, "YES", "NO")

# Evaluate performance
confusionMatrix(as.factor(ensemble_pred_v), as.factor(valid$CATs))


####### testSET2
## this creates new threshold based on test data.. 
rf_probs <- predict(fit_rf, newdata = testSET2, type = "prob")[, "YES"]
gbm_probs <- predict(fit_GBM, newdata = testSET2, type = "prob")[, "YES"]

ensemble_probs <- (rf_probs + gbm_probs) / 2

#roc_obj <- roc(response = testSET2$CATs, predictor = ensemble_probs, levels = c("NO", "YES"))

# Step 2: Find optimal threshold using Youden's J (maximizes sensitivity + specificity - 1)
#coords_result <- coords(roc_obj, x = "best", best.method = "youden", transpose = FALSE)

# View the threshold
#coords_result["threshold"]
#0.3911


# Final prediction
## Use the train or valid threshold here.   
ensemble_pred <- ifelse(ensemble_probs > 0.514, "YES", "NO")
#testSET2$CATs<-ifelse(testSET2$CATs>0.025,"YES","NO")

# Evaluate performance
confusionMatrix(as.factor(ensemble_pred), as.factor(testSET2$CATs))

###ensembled ROCs
# Create binary outcome and predicted probabilities
binary_actuals <- ifelse(testSET2$CATs == "YES", 1, 0)
predicted_probs <- as.numeric(ensemble_probs)

# Step 1: Compute ROC with 95% CI
roc_ensembled <- roc(binary_actuals, predicted_probs, ci = TRUE)
auc_ci <- ci.auc(roc_ensembled)

# Step 2: Get optimal threshold and classification
opt_threshold <- coords(roc_ensembled, "best", ret = "threshold")
predicted_class <- ifelse(predicted_probs >= opt_threshold[1,], 1, 0)

# Step 3: Confusion matrix and performance metrics
conf_mat <- table(Predicted = predicted_class, Actual = binary_actuals)
TN <- conf_mat["0", "0"]
TP <- conf_mat["1", "1"]
FN <- conf_mat["0", "1"]
FP <- conf_mat["1", "0"]

accuracy <- (TP + TN) / sum(conf_mat)
sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)

# Step 4: Plot ROC curve with annotations
tiff("PLOTS/ROC_ensembled_SET2_probCurve.tiff", width = 3, height = 3, units = "in", res = 500)  
plot(roc_ensembled, 
     col = "black", 
     lwd = 2, 
     legacy.axes = TRUE,
     xlab = "1 - Specificity", 
     ylab = "Sensitivity",
     main = paste0("ensembled (rf+gbm) Prediction ROC Curve"))

#abline(a = 0, b = 1, lty = 2, col = "gray")
# Step 5: Add text to plot
text(0.5, 0.84, labels = paste0("AUC = 0.9352"),cex=0.8)
dev.off()  

################################################################################
### ENSEMBLE_Caret method.
fitControl <- trainControl(
  method = 'cv',                   # k-fold cross validation
  number = 10,                   # number of folds
  savePredictions = 'final',    # saves predictions for optimal tuning parameter
  classProbs = TRUE, verboseIter = FALSE, # should class probabilities returned
  summaryFunction = twoClassSummary)  # results summary function


## Multiple using caretEnsembl.
# recall control parameters from above
trainControl <- trainControl(method='repeatedcv', 
                             number=10, 
                             repeats=10,
                             savePredictions='final', 
                             classProbs=TRUE,
                             summaryFunction = twoClassSummary)

algorithmList <- c('rf','gbm')

# Then run
set.seed(1234)
models <- caretList(CATs ~ ., 
                    data=train, 
                    trControl=trainControl,
                    tuneLength=10,
                    methodList=algorithmList) 
results <- resamples(models)
summary(models)
summary(results)

# plot results as a correlation
xyplot(resamples(models))
data <- modelCor(resamples(models))
data

# creating correlation matrix
library (corrplot)
corrplot(data, method = 'color', order = 'alphabet')
# prefer this one!

# Box plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(results, scales=scales)

## ensembl model of multiple
greedy_ensemble <- caretEnsemble::caretEnsemble(models)
print(summary(greedy_ensemble))

predicted_train=predict(greedy_ensemble, newdata=train)
predicted_train$pCAT <- ifelse(predicted_train$YES > 0.5, "YES", "NO") 
ens_train <- confusionMatrix(reference = as.factor(train$CATs), data = as.factor(predicted_train$pCAT), mode='everything', positive='YES')
sink("RESULTS/ENSEMBLE_train.txt")
print(ens_train)
sink()

predicted_valid=predict(greedy_ensemble, newdata=valid)
predicted_valid$pCAT <- ifelse(predicted_valid$YES > 0.5, "YES", "NO") 
ens_valid <- confusionMatrix(reference = as.factor(valid$CATs), data = as.factor(predicted_valid$pCAT), mode='everything', positive='YES')
sink("RESULTS/ENSEMBLE_valid.txt")
print(ens_valid)
sink()

# Predict on testData
predicted_test <- predict(greedy_ensemble, newdata=testSET2)
predicted_test$pCAT <- ifelse(predicted_test$YES > 0.5, "YES", "NO") 
ens_test <- confusionMatrix(reference = as.factor(testSET2$CATs), data = as.factor(predicted_test$pCAT), mode='everything', positive='YES') 

#get pvalue for ROC-AUC
testSET2$OBS <-  ifelse(testSET2$CATs== "YES", 1, 0)
roc_pvalue <- roc.area(testSET2$OBS,as.numeric(predicted_test$YES))

sink("RESULTS/ENSEMBLE_test.txt")
print(ens_test)
cat("ROC-AUC analysis\n")
print(roc_pvalue)
sink()

roc(testSET2$CATs, predicted_test$YES, ci = TRUE)

## tunning
predicted_train=predict(greedy_ensemble, newdata=train)
roc_ensembled2 <- roc(train$CATs, predicted_train$YES, ci = TRUE)

best_threshold <- coords(roc_ensembled2, "best", ret="threshold")
print(best_threshold)# 0.5577. not as good!


testSET2$OBS <-  ifelse(testSET2$CATs== "YES", 1, 0)
roc_pvalue <- roc.area(testSET2$OBS,as.numeric(predicted_test$YES))


#################end of fig 4C-E #####################################
########## Plot Fig 3B GO of significantly correlated genes ##########








###########################################################
######### 4)GO analysis####################################
go_enrichment <- enrichGO(gene          = corgenes2$genes,  # Your gene list drop CAT variable
                          keyType       = "SYMBOL",   # Key type is gene symbols
                          OrgDb         = org.Hs.eg.db,  # OrgDb for humans
                          ont           = "ALL",      # Enrichment for all GO categories
                          pAdjustMethod = "BH",       # Benjamini-Hochberg for p-value correction
                          pvalueCutoff  = 0.05,       # Adjust p-value cutoff
                          qvalueCutoff  = 0.05,       # q-value cutoff
                          readable      = TRUE)       # Convert gene IDs to gene symbols in results

# Convert GO enrichment results to a data frame
go_df <- as.data.frame(go_enrichment)
dim(go_df)

# Calculate fold enrichment using GeneRatio and BgRatio
go_df$GeneRatio_num <- sapply(go_df$GeneRatio, function(x) eval(parse(text=x)))  # Convert GeneRatio to numeric
go_df$BgRatio_num <- sapply(go_df$BgRatio, function(x) eval(parse(text=x)))  # Convert BgRatio to numeric
go_df$FoldEnrichment <- go_df$GeneRatio_num / go_df$BgRatio_num  # Calculate fold enrichment
go_df$GO_type<-go_df$ONTOLOGY
go_df_sorted <- go_df[order(go_df$FoldEnrichment),  ]  # Sort in descending order

# Step 4: Extract the top 18 GO pathways based on FoldEnrichment
top_20_go_pathways <- head(go_df_sorted, 20)

# Step 5: Print the top 18 pathways
print(top_20_go_pathways)# 

write.csv(top_20_go_pathways,"PLOTS/GO_RF_44genes.csv")
top_20_go_pathways<-read.csv("PLOTS/GO_RF_44genes.csv")
# Dot plot for GO enrichment results
# Create the dot plot
# Plot Fold Enrichment with GO types
dotplot_go<-ggplot(top_20_go_pathways, aes(x = reorder(Description, FoldEnrichment), y = FoldEnrichment, fill = -log10(p.adjust))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "gray", high = "black", name = "-log10(adj.p-value") +  # Color scale for p.adjust
  geom_text(aes(label = GO_type), hjust = -0.2, size = 2) +  # Add GO type next to bars
  labs(title = "GO Pathways by Fold Enrichment 44 genes",
       x = "GO Pathways",
       y = "Fold Enrichment",
       fill = "-log10(adj.p-value)") +
  theme_minimal() +
  theme(text = element_text(size = 8))


tiff(file = "PLOTS/dotplot_goSET1.tiff", width = 8, height = 3, units = "in", res = 500)
print(dotplot_go)
dev.off()

######### plot most important genes  #####
###########plot most important variables #################
##########################################################
varimp_rf<-varImp(fit_rf)
corelation<-corgenes2$correlation[corgenes2$genes %in% genes]

varGene<-varimp_rf[["importance"]]
names(varGene)
impgenes<-data.frame(Genes=rownames(varGene),Importance_score=varGene$Overall,corelation_ptau.Abeta42=corelation)


plot_imp <- ggplot(impgenes, aes(x = reorder(Genes, Importance_score), y = Importance_score)) + 
  geom_segment(aes(x = Genes, xend = Genes, y = 0, yend = Importance_score), color = "grey") +
  geom_point(
    aes(fill = corelation_ptau.Abeta42),  # <-- map fill to correlation
    size = 1.5, shape = 21, color = "black", stroke = 0.2
  ) +
  scale_fill_gradient(
    low = "lightblue", high = "blue", 
    name = "correlation\npTau/Aβ42", 
    guide = guide_colorbar(barheight = 2, barwidth = 0.2)
  ) +
  coord_flip() + 
  theme_minimal() + 
  labs(x = "Gene Symbols",
       y = "Importance Score of Random Forest Model") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.2, "cm"),
    
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    axis.title.y = element_text(size = 6),
    
    panel.border = element_rect(fill = NA, color = "black", size = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank()
  )

tiff(file = "PLOTS/plot_RandomForestimpVariable.tiff", width = 3, height = 3, units = "in", res = 500)
print(plot_imp)
dev.off()
#########end of most important variable  plot ################

########## corcordance plots fig 5A-D #########
# Get indices of Con patients

SET2_preds<-data.frame(preds=as.factor(output_test_SET2_rf),actual=as.factor(SET2_pheno$CAT),race=SET2_pheno$RACE,age=SET2_pheno$AGE,sex=SET2_pheno$SEX, DX=SET2_pheno$DX)

# Create a new column indicating whether the prediction was correct
SET2_preds$correct <- ifelse(SET2_preds$preds == SET2_preds$actual, "Correct", "Wrong")

SET2_preds$DX
## SUBSET data for MCI and Con patients
subset_data <- SET2_preds %>%
  filter(DX %in% c("Con", "MCI"))

concordance_summary <- subset_data %>%
  group_by(DX, correct) %>%
  summarise(count = n(), .groups = "drop")

tiff("PLOTS/concordance_ClinicalProportion.tiff", width = 3, height = 2, units = "in", res = 500)
#Plot using ggplot with lightgray & black colors
ggplot(concordance_summary, aes(x = DX, y = count, fill = correct)) +
  geom_bar(stat = "identity", position = "fill") +  # Proportional stacked bar
  scale_y_continuous(labels = scales::percent) +  # Convert y-axis to percentage
  scale_fill_manual(values = c("Correct" = "lightgray", "Wrong" = "black")) +  # Set colors
  labs(
    x = "Clinical",
    y = "Proportion",
    fill = "Prediction") +
  theme_minimal() +
  theme(text = element_text(size = 6))
dev.off()

concordance_summary
concordance_table <- matrix(c(43, 7, 24, 2), nrow = 2, byrow = TRUE,
                            dimnames = list(c("Black", "White"), c("Correct", "Wrong")))
# Perform Fisher's Exact Test
fisher_result_clinical <- fisher.test(concordance_table)
# Print the results
print(fisher_result_clinical)
######

# Fisher's Exact Test for Count Data
# 
# data:  concordance_table
# p-value = 0.71
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.04861343 3.00907855
# sample estimates:
# odds ratio 
#  0.5160207 

# Subset data for White and Black race
subset_data <- SET2_preds %>%
  filter(race %in% c("White", "Black"))
concordance_summary <- subset_data %>%
  group_by(race, correct) %>%
  summarise(count = n(), .groups = "drop")

tiff("PLOTS/concordance_RaceProportion.tiff", width = 3, height = 2, units = "in", res = 500)
#Plot using ggplot with lightgray & black colors
ggplot(concordance_summary, aes(x = race, y = count, fill = correct)) +
  geom_bar(stat = "identity", position = "fill") +  # Proportional stacked bar
  scale_y_continuous(labels = scales::percent) +  # Convert y-axis to percentage
  scale_fill_manual(values = c("Correct" = "lightgray", "Wrong" = "black")) +  # Set colors
  labs(
    x = "RACE",
    y = "Proportion",
    fill = "Prediction") +
  theme_minimal() +
  theme(text = element_text(size = 6))
dev.off()

concordance_summary
concordance_table <- matrix(c(27, 2, 37, 7), nrow = 2, byrow = TRUE,
                            dimnames = list(c("Black", "White"), c("Correct", "Wrong")))
# Perform Fisher's Exact Test
fisher_result_race <- fisher.test(concordance_table)
# Print the results
print(fisher_result_race)

# 
# Fisher's Exact Test for Count Data
# 
# data:  concordance_table
# p-value = 0.3026
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.4339133 26.7816374
# sample estimates:
# odds ratio 
#   2.524499 
#   

subset_data<-SET2_preds %>%
  filter(sex %in% c("Female", "Male"))

concordance_summary <- subset_data %>%
  group_by(sex, correct) %>%
  summarise(count = n(), .groups = "drop")

concordance_summary
concordance_table <- matrix(c(51, 4, 16, 5), nrow = 2, byrow = TRUE,
                            dimnames = list(c("Female", "Male"), c("Correct", "Wrong")))

# Perform Fisher's Exact Test
fisher_result_sex <- fisher.test(concordance_table)
# Print the results
print(fisher_result_sex)

# Fisher's Exact Test for Count Data
# 
# data:  concordance_table
# p-value = 0.1046
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.7414649 22.1961629
# sample estimates:
# odds ratio 
#    3.89728 
tiff("PLOTS/concordance_SEXProportion.tiff", width = 3, height = 2, units = "in", res = 500)
# Plot using ggplot with lightgray & black colors
ggplot(concordance_summary, aes(x = sex, y = count, fill = correct)) +
  geom_bar(stat = "identity", position = "fill") +  # Proportional stacked bar
  scale_y_continuous(labels = scales::percent) +  # Convert y-axis to percentage
  scale_fill_manual(values = c("Correct" = "lightgray", "Wrong" = "black")) +  # Set colors
  labs(
    x = "SEX",
    y = "Proportion",
    fill = "Prediction") +
  theme_minimal() +
  theme(text = element_text(size = 6))
dev.off()

###AGE######

SET2_preds$cage<-ifelse(SET2_preds$age<67.53, "Low_age", "High_age")

subset_data <- SET2_preds %>%
  filter(cage %in% c("High_age", "Low_age"))
concordance_summary <- subset_data %>%
  group_by(cage, correct) %>%
  summarise(count = n(), .groups = "drop")

concordance_summary
concordance_table <- matrix(c(32, 6, 35, 3), nrow = 2, byrow = TRUE,
                            dimnames = list(c("Female", "Male"), c("Correct", "Wrong")))

# Perform Fisher's Exact Test
fisher_result_age <- fisher.test(concordance_table)
# Print the results
print(fisher_result_age)

# Fisher's Exact Test for Count Data
# 
# data:  concordance_table
# p-value = 0.4799
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.06899139 2.37921713
# sample estimates:
# odds ratio 
#  0.4617752 
#  
# Plot using ggplot with lightgray & black colors
tiff("PLOTS/concordance_AGEProportion.tiff", width = 3, height = 2, units = "in", res = 500)
ggplot(concordance_summary, aes(x = cage, y = count, fill = correct)) +
  geom_bar(stat = "identity", position = "fill") +  # Proportional stacked bar
  scale_y_continuous(labels = scales::percent) +  # Convert y-axis to percentage
  scale_fill_manual(values = c("Correct" = "lightgray", "Wrong" = "black")) +  # Set colors
  labs(
    x = "Age Category",
    y = "Proportion",
    fill = "Prediction") +
  theme_minimal() +
  theme(text = element_text(size = 6))
dev.off()

##########end of concordance plots 5 A-D ##############




