---
title: "Predict cancer subtypes from protein level signatures with random forest model "
author: "Malina Doynova"
date: "5/02/2024"
output: github_document
  # html_notebook:
  #    theme: simplex # Change the theme to 'darkly'
  #    highlight: tango # Change the highlight style
  #    toc: yes # Include table of content
  #    toc_float: false # If to be floating or not 
  #   
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# 1.  Install and/or load the required packages

```{r}

# install.packages("BiocManager")
# BiocManager::install(c("phyloseq", "vegan", "microbiome", "Maaslin2"))
# BiocManager::install(c("microbiomeMarker"))
# BiocManager::install(c("ggplot2", "tidyverse", "dplyr", "data.table"))
# install.packages("https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-14.tar.gz", repos = NULL, type = "source")

library(mixOmics)
library(igraph)
library(randomForest)
library(caret)
library(pROC)
library(pheatmap)

```

# 2. Set the working directory

```{r}
setwd("C://Users//Malina//Desktop//RF_BC_subtypes_predict/")
```

# 3. Load the example dataset

The example dataset is from the package mixOmics. Rownames are patients and column are the proteins (sometimes with post-translational modification like "4E-BP1_pT37"). The data is normalized. 
We will call the variable mRNA, but will know the current analyses are done for normalized protein counts.
We also know what type of cancer the patients were diagnosed with, which is contained in the Y variable and will add it to the dataframe as an additional column, so we will have our outcome.  


```{r, echo = TRUE}
data(breast.TCGA)

## Explore the dataset
mRNA = breast.TCGA$data.train$protein

dim(mRNA)

head(mRNA[,1:20])

## Get Breast Cancer subtype column

Y <- breast.TCGA$data.train$subtype
summary(Y)

## Merge the protein level data with the outcome column

gene_expr_data= data.frame(cbind((mRNA),data.frame(Y)))

head(gene_expr_data[,1:20])
gene_expr_data$Y = as.factor(gene_expr_data$Y)
dim(gene_expr_data)

## Assing the dataframe with only 40 features to gene_expr_data
# gene_expr_data = mRNA_imp
# 
# dim(gene_expr_data)
# 
# head(gene_expr_data[,1:20])
# 
# levels(gene_expr_data$Y)

```


# 4. Split data into training and testing sets

```{r, echo = TRUE}

set.seed(123) # for reproducibility
train_index <- sample(1:nrow(gene_expr_data), 0.7 * nrow(gene_expr_data)) # 70% for training
train_data <- gene_expr_data[train_index, ]

test_data <- gene_expr_data[-train_index, ]

```

# 5. Train Random Forest model

```{r, echo = TRUE}
## include 10x cross-validation
ctrl <- caret::trainControl(method = "cv", number = 10)

rf_model <- randomForest(Y ~ ., data = train_data, trControl = ctrl, ntree = 500)
```

# 6. Feature Selection
```{r, echo = TRUE}
importance <- importance(rf_model)
importance = data.frame(importance)
head(importance)
importance$genes = rownames(importance)
top_features <- importance%>%dplyr::arrange(desc(MeanDecreaseGini))%>%head(40) # Select top number features
head(top_features, n=40)

```

# 7. Evaluate Model Performance and Print results

```{r, echo = TRUE}

# Model evaluation

predictions1 <- predict(rf_model, newdata = test_data)


# Compute confusion matrix
conf_matrix <- table(test_data$Y, predictions1)

# Calculate precision
precision <- conf_matrix[2, 2] / sum(conf_matrix[, 2])

# Calculate recall (sensitivity)
recall <- conf_matrix[2, 2] / sum(conf_matrix[2, ])

# Calculate specificity
specificity <- conf_matrix[1, 1] / sum(conf_matrix[1, ])

# Calculate F1 score
f1_score <- 2 * (precision * recall) / (precision + recall)

# Calculate accuracy

accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)

# Print the metrics

print(conf_matrix)
print(paste("Accuracy:", accuracy))
print(paste("Precision:", precision))
print(paste("Recall (Sensitivity):", recall))
print(paste("Specificity:", specificity))
print(paste("F1 Score:", f1_score))
```
Precision:
Precision measures the accuracy of positive predictions made by the model. It is calculated as the number of true positive predictions divided by the total number of positive predictions made by the model. A high precision indicates that the model is good at correctly identifying positive cases without incorrectly classifying negative cases as positive.

Recall (Sensitivity):
Recall, also known as sensitivity or true positive rate, measures the ability of the model to correctly identify all positive instances. It is calculated as the number of true positive predictions divided by the total number of actual positive instances in the data. High recall indicates that the model is good at capturing positive instances without missing many of them.

Specificity:
Specificity measures the ability of the model to correctly identify all negative instances. It is calculated as the number of true negative predictions divided by the total number of actual negative instances in the data. High specificity indicates that the model is good at correctly classifying negative instances without misclassifying positive instances as negative.

F1 Score:
The F1 score is the harmonic mean of precision and recall. It provides a single metric that balances both precision and recall. It is calculated as 2 times the product of precision and recall divided by the sum of precision and recall. The F1 score ranges from 0 to 1, where a higher value indicates better overall performance in terms of both precision and recall.

Accuracy:
Accuracy measures the overall correctness of the model's predictions across all classes. It is calculated as the number of correct predictions (both true positives and true negatives) divided by the total number of predictions made by the model. High accuracy indicates that the model is making correct predictions for the majority of instances.

# 8. Plot top important features averaged across the subtypes of cancer 

```{r, echo = TRUE}
cluster_genes = top_features$genes
length(cluster_genes) 

# only select the proteins that are most important features for the model 
mRNA_imp = gene_expr_data[,c(cluster_genes,"Y")]

# Calculate average expression of each gene for each cancer type

average_expr <- aggregate(. ~ Y, data = mRNA_imp, FUN = mean)

# Transpose the dataframe for better visualization
average_expr <- data.frame(t(average_expr))
colnames(average_expr) = c("Basal", "Her2", "LumA")
# Remove the 'Y' row (as it's not a gene)
average_expr <- average_expr[!rownames(average_expr) %in% "Y", ]

average_expr = data.frame(average_expr)

# Apply as.numeric to each column of the dataframe
average_expr[] <- lapply(average_expr, as.numeric)

head(average_expr)

pheatmap(t(average_expr), fontsize_col = 5,fontsize_row = 10,
         height = 18, width = 16, cluster_cols = TRUE, scale = "column")
```

From the pheatmap it is obvious there is great separation between the three cancer subtypes. Her2 pos cancer has high levels of Her2 proteins expressed. LumA has high protein levels of hormone receptors such as AR, ER, PR and also GATA3 protein levels. The Basal subtype is represented rather by proteins involved in the progression of cell cycle and those far is the most difficult to treat. 

# 9. Plot ROC and calculate AUC for each cancer subtype  
```{r, echo = TRUE}

# predict class and then attach test class
predictions <- as.data.frame(predict(rf_model, newdata = test_data, type = "prob"))
predictions$predict <- names(predictions)[1:3][apply(predictions[,1:3], 1, which.max)]
predictions$observed <- test_data$Y
head(predictions)

#  ROC curve, Basal vs non Basal
roc.Basal <- roc(ifelse(predictions$observed=="Basal", "Basal", "non-Basal"), as.numeric(predictions$Basal))
plot(roc.Basal, col = "gray60")

# others
roc.Her2 <- roc(ifelse(predictions$observed=="Her2", "Her2", "non-Her2"), as.numeric(predictions$Basal))
roc.LumA <- roc(ifelse(predictions$observed=="LumA", "LumA", "non-LumA"), as.numeric(predictions$Basal))
lines(roc.Her2, col = "blue")
lines(roc.LumA, col = "red")

roc.Basal$auc
roc.Her2$auc
roc.LumA$auc

```
