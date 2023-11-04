
# Load the libraries and data

gc() # clear memory
library(tidyverse)
library(caret)
library(randomForest)
library(rfUtilities)
getwd()

# Load the sampled presence/background data

Pres <- read_csv('data/tabular/mod/training/sampled/presence_samples.csv') %>%
  mutate(label = 1)
Bg <- read_csv('data/tabular/mod/training/sampled/background_samples.csv') %>%
  mutate(label = 0)


################################
# Set up the parallel processing

library(doSNOW)
library(doParallel)

# Set up the parallel compute
cl = makeCluster(parallel::detectCores()-1, type = "SOCK")
registerDoParallel(cl)


#############################################
# Merge the presence/background learning data

pbl <- bind_rows(Pres,Bg) %>%
  mutate(PresAbs = as.factor(label)) %>%
  dplyr::select(-c(mask,.geo,`system:index`,label,grid_id,latitude,longitude,n_sample))
glimpse(pbl)

# Clear up the memory
rm(Pres,Bg)

# Grab the dependent/independent variables
y <- pbl$PresAbs
X <- pbl %>% 
  dplyr::select(-PresAbs) %>% 
  na.omit()

# Check the % of the positive class, does it meet the 1/3 rule for sample balance?
(dim(pbl[pbl$PresAbs == 1, ])[1] / dim(pbl)[1]) * 100


###########################################################################
# Test for multicolinearity in the PBL data using the 'rfUtilities' package

set.seed(1234)

# Simple test
(mc <- rfUtilities::multi.collinear(X, perm=FALSE, p=0.05))

# Permutated test with leave out
mvm = rfUtilities::multi.collinear(
  X, perm=TRUE, leave.out=TRUE, n=101, p=0.0000001
)

# Check where frequency > 0
(rm.vars <- mvm[mvm$frequency > 0,]$variables)

# Grab a data frame with multicollinear variables removed
df <- pbl[,-which(names(pbl) %in% mvm[mvm$frequency > 0,]$variables)]

# Isolate the variables
y <- df$PresAbs
X <- df %>% dplyr::select(-PresAbs)

# Clean up
rm(mvm)


###########################
# Parameterize the RF model

# Tune the mtry parameter using randomForest
bmtry <- tuneRF(X, y, mtryStart=5, stepFactor=1.5, improve=1e-5, ntreeTry=501)
print(bmtry)

# Grab the best mtry
mtry <- data.frame(bmtry)
mtry <- mtry[which.min(mtry$OOBError),]$mtry

rm(bmtry)

#######################################################################
# Use model selection to find the most parsimonious model (rfUtilities)

rf.model <- rfUtilities::rf.modelSel(
  x=X, 
  y=y, 
  imp.scale="mir", 
  ntree=1001,
  seed=1234
)

plot(rf.model) # plot the trees

# Grab the optimal features
sel.vars <- rf.model$selvars


########################################################
# Create train/test split accounting for class imbalance

ind <- sample(2, nrow(df), replace = TRUE, prob = c(0.6, 0.4))
train <- df[ind==1,]
test <- df[ind==2,]

y <- train$PresAbs
X <- train %>% dplyr::select(-PresAbs)

# Free up space
rm(
  rf.model,df,pbl,ind
)
gc()

################################


# Now fit the random forest
rf.fit <- randomForest(
  y=y, 
  x=X[,sel.vars], # optimized features 'sel.vars'
  ntree=1001, 
  mtry=mtry, # from tuneRF
  importance=TRUE,
  proximity=TRUE,
  do.trace=50
)

# Print/plot the summaries
print(rf.fit)
plot(rf.fit)
hist(
  treesize(rf.fit),                       
  main = "Number of Nodes for the Trees",
  col = "grey")


################################

# Accuracy assessment

# Prediction & Confusion Matrix - train data
ptr <- predict(rf.fit, train[,sel.vars], type="response")
confusionMatrix(ptr, as.factor(train$PresAbs))

# Prediction & Confusion Matrix - test data
pts <- predict(rf.fit, test[,sel.vars], type="response")
confusionMatrix(pts, as.factor(test$PresAbs))

# Prepare the variable importance plot
varimp <- data.frame(importance(rf.fit, scale=TRUE, type=1)) 
varimp <- varimp %>%
  mutate(INDEX = rownames(varimp)) %>%
  arrange(desc(MeanDecreaseAccuracy))

# Plot mean decreased accuracy
impPlot <- varimp %>%
  top_n(50, MeanDecreaseAccuracy) %>%
  ggplot(aes(x = reorder(INDEX, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_bar(aes(fill=MeanDecreaseAccuracy),
           width=0.5, stat="identity", position = position_dodge(width=0.8)) +
  scale_fill_gradientn(colors = viridis_pal(begin=0.2, end=0.85, option="rocket")(3)) +
  coord_flip() +
  labs(title = "Variable Importance",
       subtitle = "Random Forests (N = 1001)",
       x= "",
       y= "Mean Decrease in Accuracy",
       caption = "") +
  theme(plot.caption = element_text(face = "italic")) +
  guides(fill="none") +
  theme_bw(14)
impPlot

# Save the model
save(rf.fit, file = "code/R/fits/rf_fit_opt.RData")
# load("code/R/fits/rf_fit_opt.RData")

# Save the feature importance plot
ggsave(impPlot, file = "figures/srme_feature_selection_rf_varImp.png",
       width=6.5, height=4.25, dpi = 300) # adjust dpi accordingly

rm(rf.fit,test,train,ptr,pts,X,varimp,impPlot)
gc() # garbage


# #########################################################
# #########################################################
# #########################################################4
# 
# # Random Forest does not account well for class imbalance
# # test the classification workflow using XGBoost
# 
# library(xgboost)
# 
# # Set up the XGBoost
# 
# # Convert outcome to logical
# train <- train %>% mutate(PresAbs = if_else(PresAbs==1,TRUE,FALSE))
# test <- test %>% mutate(PresAbs = if_else(PresAbs==1,TRUE,FALSE))
# 
# # Get into a D Matrix
# 
# #Train
# train.data <- data.matrix(train[, -which(names(train)=="PresAbs")])
# train.label <- train$PresAbs
# dtrain <- xgb.DMatrix(data = train.data, label = train.label)
# #Test
# test.data <- data.matrix(test[, -which(names(test)=="PresAbs")])
# test.label <- test$PresAbs
# dtest <- xgb.DMatrix(data = test.data, label = test.label)
# 
# # Set up the model and run it
# neg.cases = sum(train.label == FALSE)
# pos.cases <- sum(train.label == TRUE)
# 
# model.tune <- xgboost(
#   data = dtrain,
#   nround = 10,
#   max_depth = 3,
#   objective = "binary:logistic",
#   scale_pos_weight = neg.cases / pos.cases
# )
# 
# # Predict on test
# pred <- predict(model.tune,dtest)
# # Grab the error
# err <- mean(as.numeric(pred>0.5) != test.label)
# print(paste("test error:",err))
# 
# # get information on how important each feature is
# importance_matrix <- xgb.importance(names(train.data), model = model.tune)
# # and plot it!
# xgb.plot.importance(importance_matrix)
# 
# 
# # stop the cluster 
# stopCluster(cl)

