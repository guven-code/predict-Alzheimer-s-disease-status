#### Alzheimer's dementia Gaussian Mixture modeling prediction #modeling of biomarkers for the most correlated genes based on #SET1 pTau/Abeta42 levels by gene expression
#SET1 is partitioned into train(0.8) and validation(0.2).
#later tested on SET2 
##AUTHORs: EGUVEN and RMELLER
##REVISED DATE: 05/20/2025
#load required libraries
library(mclust)
#####LOAD THE PACKAGES
library(colorRamp2)
library("limma")
library("statmod")
library("edgeR")
library("sva")
library("caret")
library("MatchIt")
library("tidyr")
library("ggplot2")
library("limma")
library("edgeR")
library("MLmetrics")
library("caret")
library(caretEnsemble)
library("dplyr")
library("ggfortify")
library(MLeval)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(enrichplot)
#library(enrichplot)
library(gridExtra)
library (corrplot)
library(ComplexHeatmap)
library(grid)
library(patchwork)  
library(glmnet)
library(boot)
library(ggpubr)
library(pROC)
library(mclust)
library(Metrics)  # for MSE
library(VennDiagram)


rm(list = ls())
gc() #free up memory and report memory usage

################################################################################
##### STEP 1:  Import data and clean up tables.           ######################
################################################################################
## recode pTau/Ab to >0.0198
#setwd("~/Desktop/DESKTOP/ALZHEIMERS/EG_Pred")
df_0 <- read.csv("gene_count_matrix.csv", header = TRUE)
##Load Datasheet (Patient Phenotype data from Emory-ADRC)
pheno<- read.table("Copy of AD_Patient_data_Set1_Set2_07092024_clean.txt", sep = "\t", header = TRUE)

### Use this to change the cutoff used to call AD status.  
#pheno$CAT <- NULL
## Schindler et al 2018 used a cutoff of 0.0198 vs 0.025 from Emory ADRC
#pheno$CAT <- ifelse(pheno$"ptau.abeta" > 0.0198, "YES", "NO")

head(df_0); dim (df_0) # 62031 x 169
head(pheno); dim (pheno) # 168 x 17

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
head(exp_counts,5); dim(exp_counts)


###MATCH THE PHENO_DATA SAMPLES WITH EXP_COUNTS
# First change . to _ 
# IDs batch 1 start with an "X" so sub that out.
# IDs batch 2 start with nAD_ so sub that out
colnames(exp_counts) <- gsub("[.]","_", colnames(exp_counts))
colnames(exp_counts) <- gsub("[-]","_", colnames(exp_counts))
colnames(exp_counts) <- gsub("X","", colnames(exp_counts))
colnames(exp_counts) <- gsub("nAD_","", colnames(exp_counts))

# we have one duplicate GUID 1.65166, so remove
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
  filename = '#14_venn_diagramm.png',
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

SET2_exp_counts <- SET2_exp_counts[!grepl(smallnuclist, rownames(SET2_exp_counts)),]
dim(SET2_exp_counts) # 18923 x 76

################################################################################
####  STEP 2: Split the data into train and remaining sets based on the target variable pTau/Abeta42
################################################################################

### this keeps the same seed and index as classification modeling!
set.seed(1230)
index <- createDataPartition(SET1_pheno$CAT, p = 0.8, list = FALSE)

# Split both feature and target datasets for training data
trainData_counts <- SET1_exp_counts[, index]
#assign the train response to ptau.abeta s
trainingData_pheno <- SET1_pheno$ptau.abeta[index]
pheno_train<- SET1_pheno[index,]
dim(pheno_train)

# Split the remaining data for validation and testing
remainingData_counts <- SET1_exp_counts[, -index]
#assign the validation response to ptau.abeta s
remainingData_pheno <- SET1_pheno$ptau.abeta[-index]
pheno_valid <- SET1_pheno[-index,]
dim(pheno_valid)

validationData_counts <- remainingData_counts#[, remaining_index]
validationData_pheno <- remainingData_pheno#[remaining_index]

# Print dimensions of splits
cat("Training data dimensions:", dim(trainData_counts), "\n")
cat("Validation data dimensions:", dim(validationData_counts), "\n")

# Step 3: Transpose the count data for downstream analysis
trainData_counts <- t(trainData_counts)
validationData_counts <- t(validationData_counts)

# Names are wrong
cpm_train<-trainData_counts
cpm_valid <-validationData_counts
cpm_test <- SET2_exp_counts

################################################################################
####. STEP 3:  CREATE SIGNIFICANT GENES DATASHEET         ######################
################################################################################

# Assuming 'cpm_train' is a matrix (genes x samples) and 
# 'SET1_pheno$ptau.abeta[index]' is a vector of outcomes
correlation_results <- apply(t(trainData_counts), 1, function(x) {
  cor.test(x, SET1_pheno$ptau.abeta[index], method = "pearson")
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
results[1:50,]
summary(results$correlation)
# subset to absolute correlation > 0.5
corgenes<-results[abs(results$correlation)>0.50 ,]
corgenes2<-corgenes[complete.cases(corgenes),]
#corgenes2$correlation
corgenes2$genes <- sub("\\.cor$", "", rownames(corgenes2))
corgenes2 #check genes
dim(corgenes2)

corelation<-data.frame(cor=corgenes2$correlation,genes=corgenes2$genes)
genes<-corelation$genes

################################################################################
## Subset train data based on 44 correlated genes
## if imported list use this

train <- cpm_train[ ,colnames(cpm_train) %in% genes]
dim(train)
# Subset validation data based on 44 correlated genes
validation <- cpm_valid[ ,colnames(cpm_valid) %in% genes]
dim(validation)
str(validation)

testSet2<-cpm_test[rownames(cpm_test) %in% genes,]
str(testSet2)

################################################################################
#### STEP 4: Center and scale the training data (do this only on the training set)
################################################################################
set.seed(1230)
preProcValues <- preProcess(train, method= c("center", "scale"))
# Apply the centering and scaling transformation to training, validation, and test sets
train1 <- predict(preProcValues, train)
validation1 <- predict(preProcValues, validation)
test_SET2<-predict(preProcValues,t(testSet2))

# Check the dimensions to make sure they are correct
cat("Training set dimensions after centering and scaling:", dim(train1), "\n")
cat("Validation set dimensions after centering and scaling:", dim(validation1), "\n")
cat("Test set dimensions after centering and scaling:", dim(test_SET2), "\n")

# Prepare the final training, validation, and test data by adding the target variable CATs: ptau.abeta
# Add target variable ptau.abeta for training data
train <- data.frame(train1, CATs = SET1_pheno$pTAU.ABETA[index])  # Add CATs to the training data

# Add target variable (CATs: ptau.abeta[-index] for validation data
valid <- data.frame(validation1, CATs = SET1_pheno$pTAU.ABETA[-index])  # Add CATs to the validation data

# Add target variable (CAT) for test data
testSET2<-data.frame(test_SET2,CATs=SET2_pheno$pTAU.ABETA)

# Check final dimensions of the training, validation, and test sets
cat("Final training set dimensions:", dim(train), "\n")
cat("Final validation set dimensions:", dim(validation), "\n")
cat("SET2 dimensions:", dim(testSET2), "\n")

###################### STEP 5) START running for G=3:7 Gaussian mixture modeling to pick the optimum G ##############################
###DECISION of the optimum G
# Initialize a list to store results for different values of G
results <- list()
# Loop over different values of G (3, 4, 5, 6, 7)
results <- list()
model_metrics <- data.frame(G = integer(), logLik = numeric(), BIC = numeric(),  ICL = numeric())

# Loop over different values of G (3, 4, 5, 6, 7)
for (G in 3:7) {
  # Fit the GMM model with G components
  gmm_model <- Mclust(train$CATs, G = G)
  
  # Get the component means
  component_means <- gmm_model$parameters$mean
  
  # Predictions on train dataset
  train_pred_comp <- predict(gmm_model, train$CATs)$classification
  train_pred_CATs <- component_means[train_pred_comp]
  
  # Predictions on validation dataset
  valid_pred_comp <- predict(gmm_model, valid$CATs)$classification
  valid_pred_CATs <- component_means[valid_pred_comp]
  
  # Predictions on test dataset (SET2)
  set2_pred_comp <- predict(gmm_model, testSET2$CATs)$classification
  set2_pred_CATs <- component_means[set2_pred_comp]
  
  # Combine results into dataframes for each G
  train_df <- data.frame(actual = train$CATs, predicted = train_pred_CATs, component = factor(train_pred_comp))
  valid_df <- data.frame(actual = valid$CATs, predicted = valid_pred_CATs, component = factor(valid_pred_comp))
  set2_df  <- data.frame(actual = testSET2$CATs, predicted = set2_pred_CATs, component = factor(set2_pred_comp))
  
  # Store results in the list
  results[[paste("G", G, sep = "_")]] <- list(
    train = train_df,
    valid = valid_df,
    set2  = set2_df,
    logLik = gmm_model$loglik,
    BIC = gmm_model$BIC[2],
    ICL = gmm_model$icl
  )
  
}

# View model comparison table
print(results[c(3,4,5,6,7)])


# Data: G, BIC, ICL, and Train R²
df <- data.frame(
  G = 3:7,
  BIC = c(258, 252, 241, 236, 222),
  ICL = c(251, 240, 225, 229, 214),
  R2_Train = c(0.836, 0.9012, 0.9476, 0.9462, 0.9646)
)

# Compute BIC/ICL ratio
df$BIC_ICL_Ratio <- df$BIC / df$ICL

library(inflection)
uik(df$BIC_ICL_Ratio,df$R2_Train)

elbow_ratio <- uik(df$BIC_ICL_Ratio, df$R2_Train)
elbow_ratio
# 1.071111
elbow_point <- df[df$BIC_ICL_Ratio==elbow_ratio, ]
elbow_point
# G BIC ICL R2_Train BIC_ICL_Ratio
# 5 241 225   0.9476      1.071111
elbow_point$G
# 5


# Plot FIG6 A
tiff("PLOTS/R2_BIC.ICLratio_GMM_decisionofG.tiff", width = 3, height = 3, units = "in", res = 500)  
ggplot(df, aes(x = BIC_ICL_Ratio, y = R2_Train)) +
  geom_point(color = "black", size = 3, shape = 1) +
  geom_line(color = "blue")+# linewidth = 0.6) +
  geom_text(aes(label = G), vjust = -1.5, size = 2.5) +
  geom_point(data = elbow_point,
             aes(x = BIC_ICL_Ratio, y = R2_Train), 
             color = "red", size = 2.5) +
  annotate("text",
           x = 1.062,
           y = 0.958,
           label = "elbow point @ G =",
           color = "red", size = 2.5, fontface = "italic") +
  labs(
    title = "BIC/ICL Ratio vs R² (Train Set)",
    x = "BIC / ICL Ratio",
    y = "R² (Train)"
  ) +
  coord_cartesian(ylim = c(min(df$R2_Train), max(df$R2_Train) + 0.025)) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.6),
    panel.grid = element_blank(),  # remove grid
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    plot.title = element_text(size = 8),
    plot.margin = margin(10, 10, 10, 10)
  )
dev.off()


################################################################### STEP 6) START Gaussian Mixture Modeling  for G= 5 ################################################################
set.seed(1230)
# 1. Fit GMM to training CATs
gmm_model <- Mclust(train$CATs, G = 5)  # G = number of Gaussian components
gmm_clusters <- gmm_model$classification  # Cluster assignments
train$gmm_clusters <- gmm_clusters  # Add this as a new feature to your data

# 2. Extract component means
component_means <- gmm_model$parameters$mean

# 3. Predict component for each point in each set
# Train
train_predicted_components <- predict(gmm_model, train$CATs)$classification
train_predicted_CATs <- component_means[train_predicted_components]

# Validation
valid_predicted_components <- predict(gmm_model, valid$CATs)$classification
valid_predicted_CATs <- component_means[valid_predicted_components]

# Test SET2
testSET2_predicted_components <- predict(gmm_model, testSET2$CATs)$classification
testSET2_predicted_CATs <- component_means[testSET2_predicted_components]

# 4. Evaluate performance
cat("Train Set:\n")
print(caret::R2(train_predicted_CATs, train$CATs))
print(RMSE(train_predicted_CATs, train$CATs))
print(MSE(train_predicted_CATs, train$CATs))

cat("\nValidation Set:\n")
print(caret::R2(valid_predicted_CATs, valid$CATs))
print(RMSE(valid_predicted_CATs, valid$CATs))
print(MSE(valid_predicted_CATs, valid$CATs))

cat("\nTest SET2:\n")
print(caret::R2(testSET2_predicted_CATs, testSET2$CATs))
print(RMSE(testSET2_predicted_CATs, testSET2$CATs))
print(MSE(testSET2_predicted_CATs, testSET2$CATs))


# --- 1. Density plot with GMM fit ---
#tiff("PLOTS/GMM_train_density.tiff", width = 4, height = 3, units = "in", res = 500)  
#density_plot <- 
gmm_density <- function(x) {
  sapply(x, function(xi) {
    sum(dnorm(xi,
              mean = gmm_model$parameters$mean,
              sd = sqrt(gmm_model$parameters$variance$sigmasq),
              log = FALSE) *
          gmm_model$parameters$pro)
  })
}

# Plot with density and GMM train fit FIG 6B
tiff("PLOTS/GMM_density.tiff", width = 5, height = 3, units = "in", res = 500)  
ggplot(train_df, aes(x = actual, fill = "train")) +  # Map fill to 'train' for legend
  geom_density(color = "black", alpha = 0.4) +  # Color is fixed, but fill is mapped
  stat_function(aes(color = "GMM Fit"), fun = gmm_density, size = 0.6) +
  scale_color_manual(values = c("GMM Fit" = "blue")) +
  scale_fill_manual(name = "actual", values = c("train" = "lightgrey"), labels = c("train")) +  # Correct fill and label for legend
  labs(title = "pTau/Aβ42 density in train with GMM Fit",
       x = "pTau/Aβ42", y = "Density") + theme_bw()+
  theme(
    axis.title.x = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    plot.title = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, size     = 0.8)
  )
dev.off()

####### get the prediction and actual results: ###############
# the train_df dataframe
train_df <- data.frame(
  actual = train$CATs,  # Actual CATs values from the training set
  predicted = train_predicted_CATs,  # Predicted CATs from the GMM model
  component = train$gmm_clusters  # The component assigned to each data point
)

# the valid_df dataframe
valid_df <- data.frame(
  actual = valid$CATs,  # Actual CATs values from the training set
  predicted = valid_predicted_CATs,  # Predicted CATs from the GMM model
  component =   valid_predicted_components)# The component assigned to each data point

# the SET2_df dataframe
SET2_df <- data.frame(
  actual = testSET2$CATs,  # Actual CATs values from the training set
  predicted = testSET2_predicted_CATs,  # Predicted CATs from the GMM model
  component = testSET2_predicted_components)   # The component assigned to each data point

####### FIGURE 6C predicted values vs actual values of test data ############## 
tiff("PLOTS/GMM_SET2.tiff", width = 5, height = 3, units = "in", res = 500)  
ggplot(SET2_df, aes(x = actual, y = predicted)) +
  geom_point(color = "black", size = 3, shape = 1) +  # Points with black color
  geom_smooth(method = "lm", color = "blue", linetype = "solid", se = FALSE, size = 1) +  # Linear regression line
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # Identity line
  geom_vline(xintercept = 0.025, linetype = "dashed", color = "red", size = 1) +  # Vertical cutoff line
  geom_hline(yintercept = 0.025, linetype = "dashed", color = "red", size = 1) +  # Horizontal cutoff line
  annotate("text", x = max(SET2_df$actual) * 0.8, y = 0.022, label = "pTau/Aβ42 > 0.025", 
           color = "red", size = 3, fontface = "bold") +  # Label for horizontal line
  annotate("text", x = 0.022, y = max(SET2_df$predicted) * 0.8, label = "pTau/Aβ42 > 0.025", 
           color = "red", size = 3, angle = 90, fontface = "bold") +  # Label for vertical line
  labs(title = "SET2: Actual vs Predicted pTau/Aβ42 (GMM)",
       x = "Actual pTau/Aβ42", y = "Predicted pTau/Aβ42") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),  # Larger, centered title
    axis.title.x = element_text(size = 10),             # Axis title size
    axis.title.y = element_text(size = 10),             # Axis title size
    axis.text.x = element_text(size = 8),               # Axis text size
    axis.text.y = element_text(size = 8),               # Axis text size
    legend.title = element_text(size = 10),             # Legend title size
    legend.text = element_text(size = 8),               # Legend text size
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Border around the plot
    panel.background = element_rect(fill = "white"),    # White panel background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )
dev.off()

####FIGURE 6D
########MAKE A ROC curve convert the predicted values into binomial "YES" and "NO"s
df<-data.frame(Actual=testSET2$CATs,Predicted_GMM=as.numeric(testSET2_predicted_CATs))
df
CATs<-ifelse(testSET2$CATs>0.025,"YES","NO")
preds<-ifelse(as.numeric(testSET2_predicted_CATs)>0.025,"YES","NO")

df_bino<-data.frame(preds,CATs)
confusionMatrix(as.factor(preds),as.factor(CATs))

output_testSET2_P  <- testSET2_predicted_CATs
str(output_testSET2_P)
preds<-ifelse(as.numeric(testSET2_predicted_CATs)>0.025,"YES","NO")

# Step 1: Convert actual labels to binary (1 = YES, 0 = NO)
binary_actuals <- ifelse(CATs == "YES", 1, 0)
predicted_probs <- as.numeric(output_testSET2_P)

# Step 2: Compute ROC
roc_GMM <- roc(binary_actuals, predicted_probs)

# Step 3: Create data frame for ggplot
roc_df <- data.frame(
  Specificity = rev(roc_GMM$specificities),
  Sensitivity = rev(roc_GMM$sensitivities)
)

####FIGURE 6D
# Step 4: Plot ROC curve with ggplot2
tiff("PLOTS/ROC_GMMSET2_probCurve.tiff", width = 3, height = 3, units = "in", res = 500)  
plot(roc_GMM, 
     col = "black", 
     lwd = 2, 
     main = paste0("GMM Prediction AUC = ",as.numeric(round(roc_GMM$auc,3)), ")"),xlab = "1 - Specificity",
     ylab = "Sensitivity",
     #xlim = c(0, 1))  # Force x-axis to go from 0 to 1
     legacy.axes = TRUE)  # Makes x-axis = 1 - specificity, y-axis = sensitivity
# Add diagonal line
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()


########STEP 8)check Robustness of GMM for Supp Table 2#########
############# GMM modeling predictions ######

train <- cpm_train[,colnames(cpm_train) %in% genes]
dim(train)
# Subset validation data based on filtered genes
validation <- cpm_valid[,colnames(cpm_valid) %in% genes]

testSet2<-cpm_test[rownames(cpm_test) %in% genes, ]
str(testSet2)
# Step 1: Center and scale the training data (do this only on the training set)
preProcValues <- preProcess(train, method = c("center", "scale"))

# Step 2: Apply the centering and scaling transformation to training, validation, and test sets
train1 <- predict(preProcValues, train)
validation1 <- predict(preProcValues, validation)

testSET2<-predict(preProcValues,t(testSet2))

######### assign CATs to the pTau/Abeta levels of train(0.8) SET1
train <- data.frame(train1, CATs = SET1_pheno$pTAU.ABETA[index])  # Add CATs to the training data

# assign CATs to the pTau/Abeta levels of validation(0.2) SET1
validationData_pheno <- SET1_pheno$pTAU.ABETA[-index]
valid <- data.frame(validation1, CATs = validationData_pheno)  # Add CATs to the validation data

#assign CATs to the pTau/Abeta levels of SET2
testSET2<-data.frame(testSET2,CATs=SET2_pheno$pTAU.ABETA)

#########(i)- 10 fold repeated-CV 

set.seed(123)
folds <- createFolds(train$CATs, k = 10)

cv_results <- data.frame(
  Fold = numeric(),
  R2 = numeric(),
  RMSE = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:10) {
  # Get training and validation indices
  valid_idx <- folds[[i]]
  train_cv <- train[-valid_idx, ]
  valid_cv <- train[valid_idx, ]
  
  # Fit GMM on training fold
  gmm_model <- Mclust(train_cv$CATs, G = 5)
  component_means <- gmm_model$parameters$mean
  
  # Predict on validation fold
  pred_valid <- predict(gmm_model, newdata = valid_cv$CATs)$z %*% component_means
  
  # Save R2 and RMSE for validation
  cv_results <- rbind(cv_results, data.frame(
    Fold = i,
    R2 = caret::R2(pred_valid, valid_cv$CATs),
    RMSE = RMSE(pred_valid, valid_cv$CATs)
  ))
}

# Summarize results
cv_summary <- cv_results %>%
  summarise(
    Mean_R2 = mean(R2),
    SD_R2 = sd(R2),
    Mean_RMSE = mean(RMSE),
    SD_RMSE = sd(RMSE)
  )

print(cv_summary)

#######################################################################
####(ii) bootstrapping 100 times
set.seed(123)
n_boot <- 100  # number of bootstrap iterations

boot_results <- data.frame(
  Dataset = character(),
  Iteration = numeric(),
  R2 = numeric(),
  RMSE = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:n_boot) {
  ## Sample with replacement from training set
  boot_idx <- sample(1:nrow(train), replace = TRUE)
  boot_train <- train[boot_idx, ]
  
  ## Fit GMM model on bootstrap sample
  gmm_model <- Mclust(boot_train$CATs, G = 5)
  component_means <- gmm_model$parameters$mean
  
  ## Predict on train, valid, testSET2
  # TRAIN
  pred_train <- predict(gmm_model, newdata = train$CATs)$z %*% component_means
  boot_results <- rbind(boot_results, data.frame(
    Dataset = "Train",
    Iteration = i,
    R2 = caret::R2(pred_train, train$CATs),
    RMSE = RMSE(pred_train, train$CATs)
  ))
  
  # VALID
  pred_valid <- predict(gmm_model, newdata = valid$CATs)$z %*% component_means
  boot_results <- rbind(boot_results, data.frame(
    Dataset = "Validation",
    Iteration = i,
    R2 = caret::R2(pred_valid, valid$CATs),
    RMSE = RMSE(pred_valid, valid$CATs)
  ))
  
  # TEST SET2
  pred_test <- predict(gmm_model, newdata = testSET2$CATs)$z %*% component_means
  boot_results <- rbind(boot_results, data.frame(
    Dataset = "SET2",
    Iteration = i,
    R2 = caret::R2(pred_test, testSET2$CATs),
    RMSE = RMSE(pred_test, testSET2$CATs)
  ))
}

## Summarize bootstrapped results
library(dplyr)
boot_summary <- boot_results %>%
  group_by(Dataset) %>%
  summarise(
    Mean_R2 = mean(R2),
    SD_R2 = sd(R2),
    Mean_RMSE = mean(RMSE),
    SD_RMSE = sd(RMSE)
  )

print(boot_summary)
