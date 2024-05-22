library(caret)
library(readr)
library(dplyr)
library(randomForest)
library(ggplot2)
library(iml)

# Impute the missing values in the explanatory variables
SwissLake_Diol.imputed_RF <- rfImpute(Diol ~., ntree= 500, data = SwissLake_Diol, iter=6)

## Variable selection: removing the most correlated variables
# Remove non-numeric variables for the correlation matrix
SwissLake_Diol.imputed_RFcor <- SwissLake_Diol.imputed_RF%>%
  dplyr::select(-c(Stratification, Diol))

# Correlation matrix
cor <- cor(SwissLake_Diol.imputed_RFcor)
cor

# Find the highly correlated variables in the dataset
findCorrelation(cor, cutoff = 0.7, names = TRUE, exact = TRUE)

# Remove the most correlated variables
DiolRF_no_cor <- SwissLake_Diol.imputed_RF%>%
  dplyr::select(-c(Elevation, Cl, Ca, TN, Max.depth))

# Run the model
RF <- randomForest(Diol ~., data = DiolRF_no_cor, ntree=2000, importance = TRUE) 
RF

#Testing the number of trees plotting the obb error evolution
oob.error.data <- data.frame(
  Trees=rep(1:nrow(RF$err.rate), times=3),
  Type=rep(c("OOB", "Yes", "No"), each=nrow(RF$err.rate)),
  Error=c(RF$err.rate[,"OOB"], 
          RF$err.rate[,"N"], 
          RF$err.rate[,"Y"]))

ggplot(data=oob.error.data, aes(x=Trees, y=Error)) +
  geom_line(aes(color=Type))

##mtry tuning
x = ncol(DiolRF_no_cor)-1
resultsmtry <- vector("list", x)
for (j in 1:x) {
  
  ## 7-fold cross-validation
  # Specify the number of folds
  k <- 7
  
  # Split the data into k random folds
  folds <- sample(cut(seq(1, nrow(DiolRF_no_cor)), breaks = k, labels = FALSE))
  
  # Initialize a list to hold the results
  results <- vector("list", k)
  
  # Perform k-fold cross-validation
  for (i in 1:k) {
    
    # Split data into training and testing sets
    data.test <- DiolRF_no_cor[folds == i, ]
    data.train <- DiolRF_no_cor[folds != i, ]
    
    # Build the model using the training set
    model <- randomForest(Diol ~ . , data = data.train, mtry = j) 
    
    # Predict on the testing set
    predictions <- predict(model, data.test)
    
    # Compute performance metrics and importance
    results[[i]]$performance <- sum(predictions == data.test$Diol) / nrow(data.test) 
  }
  
  resultsmtry[[j]]$mean <- mean(unlist(results))
  resultsmtry[[j]]$SE <- sd(unlist(results))/sqrt(length(unlist(results)))
}
resultsmtry
table <- unlist(resultsmtry)
table
best_mtry=3

## Run the model with tuned parameters
RFtune <- randomForest(Diol ~., data = DiolRF_no_cor, mtry=best_mtry, ntree=4000, importance = TRUE) 
RFtune
# Extract importance
importance(RFtune)
varImpPlot(RFtune)

## 7-fold cross-validation
# Specify the number of folds
k <- 7

# Split the data into k random folds
folds <- sample(cut(seq(1, nrow(DiolRF_no_cor)), breaks = k, labels = FALSE))

# Initialize a list to hold the results
results <- vector("list", k)

# Perform k-fold cross-validation
for (i in 1:k) {
  
  # Split data into training and testing sets
  data.test <- DiolRF_no_cor[folds == i, ]
  data.train <- DiolRF_no_cor[folds != i, ]
  
  # Build the model using the training set
  model <- randomForest(Diol ~ . , data = data.train, mtry=best_mtry,importance = TRUE) 
  
  # Predict on the testing set
  predictions <- predict(model, data.test)
  
  # Compute performance metrics and importance
  results[[i]]$importance <- importance(model)
  results[[i]]$performance <- sum(predictions == data.test$Diol) / nrow(data.test) 
  
}
results

# Accumulated local effects plots
MODGlobal <- Predictor$new(RFtune, data = DiolRF_no_cor, type = "prob", class = 2)
MAATplotG <- plot(FeatureEffect$new(MODGlobal, feature = "MAAT"))
plot(MAATplotG) +
  theme(axis.title=element_blank()) +
  theme(axis.text= element_text(size = 12)) +
  geom_hline(yintercept=0, linetype="dashed", color = "deepskyblue4", linewidth = 0.8)