
# Load the libraries and data

library(tidyverse)
library(sf)
library(caret)
library(randomForest)
library(rfUtilities)
library(ggcorrplot)

getwd()


#############################################
# Prep the sampled presence/background data #
#############################################

# Pres.path <- 'data/spatial/mod/training/sampled/presence/pi_points_srme_m500_sampled_.csv'
# Abs.path <- 'data/spatial/mod/training/sampled/background/bg_points_evt_sbcls_sampled_.csv'

Pres.path <- 'data/tabular/mod/training/presence/southern_rockies_pres_full_sampled.csv'
Abs.path <- 'data/tabular/mod/training/background/southern_rockies_bg_prop_sampled.csv'

# Load the presence data
Pres <- read_csv(Pres.path) %>%
  as_tibble() %>%
  mutate(label = 1) %>%
  na.omit()

# Load the absence/background data
Abs <- read_csv(Abs.path) %>%
  as_tibble() %>%
  mutate(label = 0) %>%
  na.omit()

# Merge and tidy the PresAbs data
PresAbs <- bind_rows(Pres,Abs) %>%
  mutate(PresAbs = as.factor(label)) %>%
  dplyr::select(-c(mask,.geo,))
  # dplyr::select(-c(`system:index`,EVT_SBCLS,SBCLS_CODE,Block_ID,fid,label,.geo))
glimpse(PresAbs)

PresAbs <- PresAbs %>% select(-c(elevation,slope))

# Clear up the memory
rm(Pres,Abs)

# Grab the dependent/independent variables
y <- PresAbs$PresAbs
X <- PresAbs %>% dplyr::select(-PresAbs) 

# Check the % of the positive class, does it meet the 1/3 rule for sample balance?
if ((dim(PresAbs[PresAbs$PresAbs == 1, ])[1] / dim(PresAbs)[1]) * 100 <= 33.0) {
  print("Lower than the 1/3 rule ...")
  print((dim(PresAbs[PresAbs$PresAbs == 1, ])[1] / dim(PresAbs)[1]) * 100)
}


#################################################################################
# Test for multicolinearity in the PresAbs data using the 'rfUtilities' package #

set.seed(1234)

# Correlation Matrix Plot

thresh <- 0.6
corr <- cor(X)
corr[abs(corr) < thresh] <- 0

ggcorrplot(corr, method = "square", 
           type = "lower", 
           lab = FALSE, 
           lab_size = 3, 
           colors = c("#6D9EC1", "white", "#E46726"), 
           title = "Correlation Matrix",
           ggtheme = theme_minimal())

# Identify the highly correlated pairs
corr.vars <- which(abs(corr) > thresh, arr.ind = TRUE)
# Removing self-correlations (diagonal elements)
(corr.vars <- corr.vars[corr.vars[,1] != corr.vars[,2], ])

# Print feature pairs
if(nrow(corr.vars) > 0) {
  for (i in 1:nrow(corr.vars)) {
    print(paste(names(X)[corr.vars[i, 1]], 
                names(X)[corr.vars[i, 2]], 
                sep=" - "))
  }
} else {
  print("No pairs with correlation above the threshold")
}

##############################
# Robust multicolinearity test

# Permutated test with leave-one-out
cl.test = rfUtilities::multi.collinear(X, perm=TRUE, leave.out=TRUE, n=999, p=0.05)
# Check where frequency > 0
(rm.vars <- cl.test[cl.test$frequency > 0,]$variables)

# Remove collinear 
df <- PresAbs[,-which(names(PresAbs) %in% rm.vars)]
df <- df %>% na.omit()

# Isolate the variables
y <- df$PresAbs
X <- df %>% dplyr::select(-PresAbs) 

# Clean up
rm(corr,corr.vars,cl.test,i,thresh)


# #############################
# # Parameterize the RF model #
# #############################
# 
# # Tune the mtry parameter using randomForest
# bmtry <- tuneRF(X, y, mtryStart=2, stepFactor=1.5, improve=1e-5, ntreeTry=999)
# print(bmtry)
# 
# # Grab the best mtry
# mtry <- data.frame(bmtry)
# (mtry <- mtry[which.min(mtry$OOBError),]$mtry)
# 
# rm(bmtry)


##########################################################
# Create train/test split accounting for class imbalance #
##########################################################

rm(y,X)

ind <- sample(2, nrow(df), replace = TRUE, prob = c(0.7, 0.3))
train <- df[ind==1,]
test <- df[ind==2,]

y_train <- train$PresAbs
y_train <- factor(y_train, levels = c(0, 1))
X_train <- train %>% dplyr::select(-PresAbs)

y_test <- test$PresAbs
y_test <- factor(y_test, levels = c(0, 1))
X_test <- test %>% dplyr::select(-PresAbs)

# Free up space
rm(df,PresAbs,ind,train,test)
gc()


###############################
# Fit the random forest model #
###############################

mtry.sqrt <- as.integer(sqrt(ncol(X_train)))

# Now fit the random forest
rf.fit <- randomForest(
  y=y_train, 
  x=X_train, # highly correlated bands removed
  ntree=1001, 
  mtry=6, # from tuneRF
  importance=TRUE
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

# Prediction & Confusion Matrix - test data
preds <- predict(rf.fit, X_test, type="response")
preds <- factor(preds, levels = levels(y_test))

confusionMatrix(preds, as.factor(y_test))

# Prepare the variable importance plot
varimp <- data.frame(importance(rf.fit, scale=TRUE, type=1)) 
varimp <- varimp %>%
  mutate(INDEX = rownames(varimp)) %>%
  arrange(desc(MeanDecreaseAccuracy))

# Plot mean decreased accuracy
impPlot <- varimp %>%
  top_n(10, MeanDecreaseAccuracy) %>%
  ggplot(aes(x = reorder(INDEX, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_bar(aes(fill=MeanDecreaseAccuracy),
           width=0.5, stat="identity", position = position_dodge(width=0.8)) +
  scale_fill_gradientn(colors = viridis::viridis_pal(begin=0.2, end=0.85, option="rocket")(3)) +
  coord_flip() +
  labs(title = "Variable Importance",
       subtitle = "Random Forests (NTree = 101)",
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
ggsave(impPlot, file = "figures/rf-tuning_best_model_feat_imps.png",
       width=6.5, height=4.25, dpi = 300) # adjust dpi accordingly

rm(rf.fit,test,train,X,varimp,impPlot)
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

