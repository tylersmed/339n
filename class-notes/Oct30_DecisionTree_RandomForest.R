library (dslabs)
str(brca)
?brca

#target variable brca$y

brca$y
set.seed(3)

#install.packages("caret")
library (caret)
library (ggplot2)
library (cowplot)

brca_df = data.frame(Target=brca$y, Feature= brca$x)


## In the previous class we used the sample function
?sample
inTrain <- sample (1:nrow(brca_df), replace=F, size= round(nrow(brca_df)*0.9))
training = brca_df[inTrain, ]
testing = brca_df[-inTrain,]
str(training)
str(testing)

## Let's use createDataPartition() function from the caret package 
## for a more robust sampling strategy.
## to maintain the ratio of target classes in both training and testing sets.
## Split data into training and testing sets
set.seed(3)
inTrain <- createDataPartition(y = brca_df$Target, p = 0.9, list = FALSE)
training <- brca_df[inTrain, ]
testing <- brca_df[-inTrain, ]

#Question-1: createDataPartition() 
#https://www.investopedia.com/terms/stratified_random_sampling.asp#:~:text=Stratified%20random%20sampling%20is%20a,as%20income%20or%20educational%20attainment.

# Create visualizations for the training data
# Histograms and boxplots by feature and target class
cowplot::theme_cowplot()
p1 <- ggplot (training, aes(x= Feature.radius_mean, fill= Target, color= Target)) + geom_histogram(binwidth=1) + theme_cowplot(12)
p2 <-ggplot (training, aes(x=Feature.perimeter_mean, fill= Target, color= Target)) + geom_histogram(binwidth=1)+ theme_cowplot(12)
p3 <-ggplot (training, aes(x=Feature.symmetry_worst, fill= Target, color= Target)) + geom_boxplot() + theme_cowplot(12)
p4 <-ggplot (training, aes(x=Feature.symmetry_se, fill= Target, color= Target)) + geom_boxplot() + theme_cowplot(12)

# Combine the plots into a grid
plot_grid(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), label_size = 12)

#install.packages("tree")
library(tree)  

#single split
decision_tree1 = tree(Target~Feature.radius_mean, data=training, control=tree.control(nrow(training), mincut=200))
summary(decision_tree1)
plot(decision_tree1)
text(decision_tree1, pretty=0)

#multi-feature tree
decision_tree3= tree(Target~ ., data=training )
summary(decision_tree3)
plot(decision_tree3)
text(decision_tree3, pretty=0)

tree.pred = predict(decision_tree3, testing, type="class")
table(tree.pred, testing$Target)

head (brca_df)
brca_df$Target
training_two = brca_df[1:360,]
testing_two = brca_df[360:562,]

dim(training)
dim(training_two)

#train a decision tree
decision_tree4= tree(Target~., data=training_two)
summary(decision_tree4)

summary(decision_tree3)
tree.pred = predict(decision_tree3, testing, type="class")
sum (tree.pred == testing$Target)/nrow(testing)

tree.pred2 = predict(decision_tree4, testing_two, type="class")
sum(tree.pred2 ==testing_two$Target)/nrow(testing_two)

#What happened? Q2

decision_tree5= tree(Target~Feature.radius_mean + Feature.symmetry_mean, data=training, 
                     control = tree.control(nrow(training), minsize = 2, mindev=0))
## minsize = 2: the smallest size of any terminal node (leaf) should be 2. 
## the algorithm will not attempt to split a node if it contains fewer than 2 cases.
#mindev = 0 sets the minimum decrease in the deviance required to attempt any split. 
## Any decrease is sufficient to make a split: 
## it does not require a specific minimum amount of improvement in model fit 
## to justify a split.
decision_tree6 = tree(Target ~Feature.radius_mean + Feature.symmetry_mean, data=training)

summary(decision_tree5)
summary(decision_tree6)

plot (decision_tree5)
plot (decision_tree6)

tree.pred3 = predict(decision_tree5, testing, type="class")
tree.pred4 = predict(decision_tree6, testing, type="class")

sum(tree.pred3 == testing$Target)/nrow(testing)
sum(tree.pred4 == testing$Target)/nrow(testing)
# the simpler tree (decision_tree6) gave better predictions on the testing set.
#Why? Q3

predict(decision_tree6, testing, type="class")
predict(decision_tree6, testing, type="vector")
# result as a category vs probabilistic
table(tree.pred4, testing$Target)

#false positive rate FPR -- FP/(FP+TN)
# true positive Rate /Sensitivity-- TP/FP+FN 
# true negative rate /Specificity -- TN/TN+FP
#precision -positive predictive value -- TP/TP + FP

#install.packages("randomForest")
library(randomForest)

rf1 = randomForest(Target ~ ., data= training, importance =T)
rf1
rf_predict = predict(rf1, testing)
sum (rf_predict == testing$Target) /nrow(testing)
#use the importance() function to view the importance measures.
## We discussed Gini index in the class
## MeanDecreaseGini reflects how much a particular variable contributes to the 
# homogeneity of the nodes and leaves in the resulting Random Forest model. 
# A higher value of MeanDecreaseGini indicates a variable is more important 
# for the model predictions.

importance(rf1) 
#importance = TRUE doesn't affect model performance
# it only calculates additional post-hoc statistics.

#visualize importance
varImpPlot(rf1)


#Q4: MeanDecreaseGini, VarImpPlot

## Using caret to tune the model
## tuning parameters such as the number of trees, the depth of trees, 
## and the number of features considered at each split.
## using the train function from the caret package with 10-fold cross-validation 
## What is cross-validation and why it might not be as relevant for RandomForest algorithm?
## Q5 : Crossvalidation
control <- trainControl(method="cv", number=10)
tuneGrid <- expand.grid(.mtry=c(1:sqrt(ncol(training)-1)))
rfTuned <- train(Target~., data=training, method="rf", trControl=control, tuneGrid=tuneGrid)
rfTuned

rfTuned_predict = predict(rfTuned, testing)
sum (rfTuned_predict == testing$Target) /nrow(testing)




## Partial Dependence Plots show the dependence 
## between the target function and a set of ‘interesting’ features
# Using the pdp package for partial dependence plots
#install.packages("pdp")
library(pdp)
partial1 <- partial(rfTuned, pred.var = "Feature.perimeter_worst", train = training)
autoplot(partial1)
dev.off()
#"y-hat," represents the predicted or estimated value of the target variable
partial2 <- partial(rfTuned, pred.var = "Feature.radius_mean", train = training)
autoplot(partial2)

# Using randomForest for clustering
rfClust <- randomForest(x=training[-ncol(training)], y=NULL, proximity=TRUE)
# Since y=NULL, there is no target variable,function performs unsupervised learning.
#  proximity=TRUE computes a proximity matrix, which measures the similarity 
# between each pair of samples in the data
# two samples are more similar if they end up in the same terminal nodes of trees.
str(rfClust)


#calculate dissimilarity from proximity scores
dissimilarity <- 1 - rfClust$proximity

#perform hierarchical clustering on the dissimilarity matrix
hc <- hclust(as.dist(dissimilarity), method="average")
plot(hc)

# Q6: What is the purpose of setting y=NULL in this context?

##Discussion: How do we choose between using Naive Bayes classification and Random Forest?
## Which factors might influence our decision 

##Random Forest is often used in bioinformatics for tasks like predicting protein-protein 
## interactions, gene selection for disease classification, or biomarker discovery 
## because of its ability to handle noisy and complex data.
## Naive Bayes is used in genomics for tasks like gene expression classification, 
## patient diagnosis based on genetic variation, or phylogenetics classification etc..
## Difference in efficiency, interpretability, data size, class imbalance
#Example: Use Naive Bayes when When dealing with high-dimensional data where features are independent and computational efficiency is critical.

