# ---- PRELIMINARIES ----
# Load necessary libraries and explore the BRCA dataset
#install.packages("dslabs")
library(dslabs)
str(brca)
?brca

# ---- DATA PREPARATION ----
# Define target variable and create a new dataframe
brca$y
set.seed(3) # Setting a seed for reproducibility
brca_df = data.frame(Target = brca$y, Feature = brca$x)

# Splitting data into training and testing sets
inTrain <- sample(1:nrow(brca_df), replace = FALSE, size = round(nrow(brca_df) * 0.9))
training <- brca_df[inTrain, ]
testing <- brca_df[-inTrain,]
str(training)
str(testing)


# ---- VISUALIZATION ----
# Load additional libraries for visualization
#install.packages("ggplot2")
#install.packages("cowplot")
library(ggplot2)
library(cowplot)

# Create visualizations for the training data
# Histograms and boxplots by feature and target class
p1 <- ggplot (training, aes(x= Feature.radius_mean, fill= Target, color= Target)) + geom_histogram(binwidth=1) + theme_cowplot(12)
p2 <-ggplot (training, aes(x=Feature.perimeter_mean, fill= Target, color= Target)) + geom_histogram(binwidth=1)+theme_cowplot(12)
p3 <-ggplot (training, aes(x=Feature.symmetry_worst, fill= Target, color= Target)) + geom_boxplot() + theme_cowplot(12)
p4 <-ggplot (training, aes(x=Feature.symmetry_se, fill= Target, color= Target)) + geom_boxplot() + theme_cowplot(12)

# Combine the plots into a grid
plot_grid(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), label_size = 12)
malignant_radius = median (training$Feature.radius_mean[training$Target == "M"])
begnin_radius = median (training$Feature.radius_mean[training$Target == "B"])

#Instapoll: Training and Visualization
## CAN WE PREDICT BASED ON TRAINING DATASET?


predictions= c()
for (example in 1:nrow(training)){
  if (training$Feature.radius_mean[example]> malignant_radius)
  {
    predictions[example]="M"
  }else {
    predictions[example]="B"
  }
}

##Evaluation 
predictions == training$Target
sum(predictions == training$Target)/nrow(training)
##80%

table(predictions)
table(training$Target)
#print the confusion matrix
table(predictions, training$Target)

naive_predictions_1 = rep("B", nrow(training))
sum(naive_predictions_1 == training$Target)/nrow(training)
##62% overall accuracy

table (naive_predictions_1, training$Target)

naive_predictions_2 = sample(c("M","B"), size=nrow(training), prob= c(196,316), replace=T)
sum(naive_predictions_2==training$Target)/nrow(training)

table(naive_predictions_2, training$Target)
# did we select the best possible threshold?
# can we leverage the other variables?
# how about playing training/test ?
#Instapoll: Naive Bayes

# ---- NAIVE BAYES CLASSIFICATION  
# Establish a baseline using naive predictions
# (Several approaches: predictions, naive_predictions_1, naive_predictions_2)
# Evaluation of the baseline classifiers
#install.packages("caret")
library (caret)
#When you train a Naive Bayesß classifier using a labeled dataset, 
#the model learns the probability of a data point belonging to each 
#category based on the feature

#install.packages("e1071")
library(e1071)  # This library contains the naiveBayes function
nb_model = naiveBayes(Target ~ ., data = training)
#naiveBayes function from the e1071 package in R can work with numerical data. 
# applies Bayes' theorem with naive independence assumptions between the features. 
#It can handle both categorical and numerical data, but it treats all variables independently given the class variable.

nb_predict = predict(nb_model, testing)
nb_accuracy = sum(nb_predict == testing$Target) / nrow(testing)  # Calculate accuracy
nb_predict  # The predictions from the Naive Bayes model

#kNN 
library(class)
train_features <- training[ , -1]  # Exclude the target column
test_features <- testing[ , -1]
train_target <- training$Target

# Use kNN to predict on the test set (k = 3 as an example)
knn_predictions <- knn(train = train_features, test = test_features, cl = train_target, k = 3)
knn_accuracy <- sum(knn_predictions == testing$Target) / length(testing$Target)
print(knn_accuracy)
table(knn_predictions, testing$Target)

# Try different values of k
k_values <- c(1, 3, 5, 7, 10)
accuracies <- c()

for (k in k_values) {
  knn_predictions <- knn(train = train_features, test = test_features, cl = train_target, k = k)
  knn_accuracy <- sum(knn_predictions == testing$Target) / length(testing$Target)
  accuracies <- c(accuracies, knn_accuracy)
  print(paste("k =", k, "Accuracy =", knn_accuracy))
}

# Plot the accuracy against k values
plot(k_values, accuracies, type = "b", xlab = "k Value", ylab = "Accuracy", main = "Accuracy vs k")

### kNN typically uses Euclidean distance as the 
# default, but you can also use other metrics like Manhattan or Minkowski distance
#  by modifying the dist argument in the kNN function.

#install.packages("caret")
library(caret)

# Define training control with kNN
train_control <- trainControl(method = "cv", number = 10)  # 10-fold cross-validation

# Train the kNN model with a different distance metric (e.g., Manhattan)
set.seed(123)
knn_model <- train(Target ~ ., data = training,
                   method = "knn",
                   trControl = train_control,
                   tuneGrid = expand.grid(k = c(1, 3, 5, 7, 10)),
                   preProcess = c("center", "scale"),
                   metric = "Accuracy",
                   tuneLength = 10)

print(knn_model)

# Use caret package for cross-validation
set.seed(123)
knn_cv_model <- train(Target ~ ., data = training,
                      method = "knn",
                      trControl = trainControl(method = "cv", number = 10),
                      tuneGrid = expand.grid(k = 3:10))

print(knn_cv_model)
plot(knn_cv_model)


#By default, kNN gives equal weight 
#to all neighbors. Let's use distance-weighted 
#voting, where closer neighbors are given more 
#influence than farther ones by specifying weight = "distance".

# Using the kNN Package for Weighted Voting
# Train a kNN model using distance weighting
#install.packages("kknn")
library(kknn)
kknn_model <- train.kknn(Target ~ ., data = training,
                         kmax = 10, distance = 2, kernel = "optimal")

kknn_predictions <- predict(kknn_model, newdata = testing)

kknn_accuracy <- sum(kknn_predictions == testing$Target) / nrow(testing)
kknn_accuracy


##Instapoll kNN questions 

# k-means clustering 
#install.packages("factoextra")
#install.packages("cluster")

library(factoextra)
library(cluster)
# Data scaling, determining optimal number of clusters, and k-means clustering
# Visualization of clustering results
brca_df = data.frame(Target=brca$y, Feature= brca$x)
str(brca_df)
brca_df <- na.omit(brca_df)
brca_scaled= scale(brca_df[,-1])

#kmeans(data, centers, nstart)
#centers: The number of clusters, denoted k.
#nstart: The number of initial configurations. 
#Because it’s possible that different initial starting clusters can lead to different results,
#it’s recommended to use several different initial configurations. 
#The k-means algorithm will find the initial configurations that lead to the smallest within-cluster variation.

# create a plot of the number of clusters vs. the total within sum of squares:
#creating a visualization of the number of clusters (k) versus
#the total within-cluster sum of squares (WSS). 
#The fviz_nbclust function is part of the factoextra package, which provides 
#functions to extract and visualize the output of 
#multivariate data analyses, including clustering. 
#This plot helps to find the optimal number of clusters using the elbow method. 

fviz_nbclust(brca_scaled, kmeans, method = "wss")
#The Total Within Sum of Square (WSS) is a measure of the internal 
#coherence of the clusters in a clustering algorithm like k-means.

# Elbow method looks at the percentage of variance explained as a function 
# of the number of clusters: one should choose a number of clusters so that 
#adding another cluster doesn’t give much better modeling of the data.
#cluster package  
gap_stat <- clusGap(brca_scaled,
                    FUN = kmeans,
                    nstart = 25,
                    K.max = 10,
                    B = 50)
#The gap statistic is a method for estimating the number of clusters in a set of data. 
#It compares the total within intra-cluster variation for different values of k with 
#their expected values under null reference distribution of the data
# The highest gap statistic suggests the optimal number of clusters

fviz_gap_stat(gap_stat)
# Let's perform k-means clustering using k of 2
set.seed(1)
#perform k-means clustering on the scaled dataset 
km <- kmeans(brca_scaled, centers = 2, nstart = 25)
km
#Let's visualize the clusters formed by k-means using the fviz_cluster function. 
fviz_cluster(km, data = brca_scaled)
brca_df[430,1]
brca_df[300,1]
str(km)
#Let's calculate the mean of each variable for each cluster. 
#use the aggregate function to apply the mean function to the brca_scaled dataset, 
#which is grouped by the cluster assignment from the k-means result.
aggregate(brca_scaled, by=list(cluster=km$cluster), mean)
km$cluster

table <- table(km$cluster, brca_df$Target)
table

#Determine the most frequent label in each cluster.
cluster_to_label_map <- apply(table, 1, function(row) {
  if (row['B'] > row['M']) {
    return('B')
  } else {
    return('M')
  }
})

# Use the mapping to convert cluster numbers to predicted labels.
predicted_labels <- factor(sapply(km$cluster, function(cluster_number) {
  cluster_to_label_map[as.character(cluster_number)]
}), levels = levels(brca_df$Target))

# cross-tabulation of cluster assignments and true labels
table <- table(km$cluster, brca_df$Target)
print(table)

# Calculate "accuracy"
accuracy <- sum(predicted_labels == brca_df$Target) / length(brca_df$Target)
print(accuracy)

#instapoll: Kmeans questions, kNN/kMeans



