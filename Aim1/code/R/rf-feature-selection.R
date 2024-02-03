
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

# Pres.path <- 'data/tabular/mod/training/presence/southern_rockies_pres_full_sampled.csv'
# Abs.path <- 'data/tabular/mod/training/background/southern_rockies_bg_prop_sampled.csv'

Pres.path <- 'data/tabular/mod/training/presence/southern_rockies_pres_full_sampled_n55.csv'
Abs.path <- 'data/tabular/mod/training/background/southern_rockies_bg_prop_sampled_n55.csv'

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
  dplyr::select(-c(label,mask,.geo,n_sample,grid_id,`system:index`,latitude,longitude))
  # dplyr::select(-c(`system:index`,EVT_SBCLS,SBCLS_CODE,Block_ID,fid,label,.geo))
glimpse(PresAbs)

# Not sure why these columns persisted ...
PresAbs <- PresAbs %>% 
  select(-c(ChlRE_autumn,ChlRE_summer,NDMI_autumn,NDMI_summer,NDRE_autumn,NDRE_summer)) %>%
  na.omit() # check for any other NA values

# # Version without elevation
# PresAbs <- PresAbs %>% select(-c(elevation,slope,aspect))

# Clear up the memory
rm(Pres,Abs)

# Grab the dependent/independent variables
y <- PresAbs$PresAbs
X <- PresAbs %>% dplyr::select(-PresAbs) 

# Check the % of the positive class, does it meet the 1/3 rule for sample balance?
if ((dim(PresAbs[PresAbs$PresAbs == 1, ])[1] / dim(PresAbs)[1]) * 100 <= 33.0) {
  paste0("Ratio of background:presence - ",
         paste0(round((dim(PresAbs[PresAbs$PresAbs == 1, ])[1] / dim(PresAbs)[1]) * 100), ":1"))
}


#################################################################################
# Test for multicolinearity in the PresAbs data using the 'rfUtilities' package #

set.seed(1234)

# Correlation Matrix Plot

thresh <- 0.8
corr <- cor(X)
corr[abs(corr) < thresh] <- 0

ggcorrplot(corr, method = "square", 
           type = "lower", 
           lab = FALSE, 
           lab_size = 3, 
           colors = c("#6D9EC1", "white", "#E46726"), 
           title = "Correlation Matrix",
           ggtheme = theme_minimal())


###################################
# Robust multicolinearity test(s) #

# Single test with no permutation
cl <- rfUtilities::multi.collinear(X, perm=FALSE, p=0.05)
for(l in cl) { # test these variables using "leave one out"
  cl.test <- X[,-which(names(X)==l)]
  print(paste("Remove variable", l, sep=": "))
  multi.collinear(cl.test, p=0.05)
}

# Permutated test with leave-one-out
cl.perm = rfUtilities::multi.collinear(
  X, perm=TRUE, leave.out=TRUE, n=1001, p=0.05
)
# Check the histogram of the frequency > 0
hist(cl.perm$frequency[cl.perm$frequency > 0])

# Check where frequency > 10% of permutations
th <- 1001*0.10 # 10% of the permutations
(rm.vars <- cl.perm[cl.perm$frequency > th,]$variables)

# Remove collinear variables
df <- PresAbs[,-which(names(PresAbs) %in% rm.vars)]
df <- df %>% na.omit()

# Isolate the variables
y <- df$PresAbs
X <- df %>% dplyr::select(-PresAbs) 

# Clean up
rm(corr,corr.vars,cl.test,thresh,Abs.path,Pres.path,l,th,cl.perm,cl)

# Check the corrplot again
thresh <- 0.8
corr <- cor(X)
corr[abs(corr) < thresh] <- 0
ggcorrplot(corr,
           type = "lower", 
           lab = FALSE, 
           lab_size = 3, 
           colors = c("#6D9EC1", "white", "#E46726"), 
           title = "Correlation Matrix",
           ggtheme = theme_minimal())

rm(corr,thresh,rm.vars,PresAbs)


###########################
# Perform model selection #
###########################

rf.model <- rf.modelSel(
  x=X, y=y,imp.scale="mir",ntree=1001,seed=42,parsimony=0.03
)

(sel.vars <- rf.model$selvars)

plot(rf.model)              # plot importance for selected variables
plot(rf.model, imp = "all") # plot importance for all variables 

# Look at the number of features for each model tested
rf.model$test

# Check one of the alternative sel.vars
(sel.vars_alt <- rf.model$parameters[[2]])

# Create a new data frame for the optimized model coefficients

X.sel <- X[,sel.vars]

rm(rf.model)


##################################
# Parameterize the mtry argument #
##################################

# Tune the mtry parameter using randomForest
bmtry <- tuneRF(
  X.sel, y, 
  mtryStart=2, 
  stepFactor=1.5, 
  improve=0.01, 
  ntreeTry=101
)
print(bmtry)

# Grab the best mtry
mtry <- data.frame(bmtry)
(mtry <- mtry[which.min(mtry$OOBError),]$mtry)

rm(bmtry)

# Or take the square root of the length of X

mtry.sqrt <- as.integer(sqrt(ncol(X.sel)))


###########################
# Create train/test split #
###########################

rm(y,X,X.sel)

ind <- sample(2, nrow(df), replace = TRUE, prob = c(0.7, 0.3))
train <- df[ind==1,]
test <- df[ind==2,]

y_train <- train$PresAbs
y_train <- factor(y_train, levels = c(1,0))
X_train <- train %>% dplyr::select(-PresAbs)
X_train <- X_train[,sel.vars]

y_test <- test$PresAbs
y_test <- factor(y_test, levels = c(1,0))
X_test <- test %>% dplyr::select(-PresAbs)
X_test <- X_test[,sel.vars]

# Free up space
rm(df,PresAbs,ind,train,test)
gc()


###############################
# Fit the random forest model #
###############################

# Now fit the random forest
rf.fit <- randomForest(
  y=y_train, 
  x=X_train, # highly correlated bands removed
  ntree=1001, 
  mtry=mtry.sqrt, # from tuneRF
  importance=TRUE
)

# Print/plot the summaries
print(rf.fit)
plot(rf.fit)
hist(
  treesize(rf.fit),                       
  main = "Number of Nodes for the Trees",
  col = "grey")


#######################
# Accuracy assessment #
#######################

# Prediction on the testing data
preds <- predict(rf.fit, X_test, type="response")
preds <- factor(preds, levels = levels(y_test))

# Print the confusion matrix
(cm <- confusionMatrix(preds, as.factor(y_test)))

# Look at the precision, recall, and F1-score
recall <- cm$byClass['Sensitivity']  # Recall is Sensitivity
specificity <- cm$byClass['Specificity']
precision <- cm$byClass['Pos Pred Value']  # Precision calculation (Positive Predictive Value
# F1 Score calculation
F1 <- 2 * (precision * recall) / (precision + recall)

# Printing the metrics
print(paste("Precision:", precision))
print(paste("Recall (Sensitivity):", recall))
print(paste("F1 Score:", F1))


#######################
# Variable Importance #
#######################

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
       subtitle = "Random Forests (NTree = 1001)",
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

rm(rf.fit,X_test,X_train,varimp,impPlot)
gc() # garbage


