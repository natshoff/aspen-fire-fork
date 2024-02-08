
# LIBRARIES

library(tidyverse)

# Laod the data (testing parts from model folds)

# testPart <- 'data/tabular/mod/results/best_model/southern_rockies_test_probs.csv'

testPart <- "data/tabular/mod/results/best_model/southern_rockies_test_probs_n55.csv"

# Sampled testing (probabilities)
test <- read_csv(testPart) %>%
  rename(gee_id = `system:index`,
         TrueLabel = label) %>%
  dplyr::select(-.geo) %>%
  mutate(probability = probability*0.001) # scale back
glimpse(test)


##################################
# Calculate the accuracy metrics #

models <- unique(test$seed)

# create a list to hold the results
i=1
ls <- list()
for(model in models){
  
  # subset to the seed #
  test.part <- test %>% filter(seed == model)
  
  # grab presence/background
  pres <- test.part %>% filter(TrueLabel==1)
  abs <- test.part %>% filter(TrueLabel==0)
  
  # get the df sizes
  pres.size <- nrow(pres)
  abs.size <- nrow(abs)
  
  # loop through cutoff values
  vec.list <- list()
  x = 1
  for(cutoff in seq(0,1,length=100)){
    
    # Create an empty vector to store results
    vec <- numeric(100)
    
    # Create the confusion matrix
    tp = pres %>% filter(probability >= cutoff) %>% nrow() 
    tpr = tp / pres.size 
    fn = pres %>% filter(probability < cutoff) %>% nrow() 
    tn = abs %>% filter(probability < cutoff) %>% nrow() 
    tnr = tn / abs.size 
    fp = abs %>% filter(probability >= cutoff) %>% nrow() 
    fpr = fp / abs.size 
    
    # Assign to vector and append
    vec <- c(model,pres.size,abs.size,cutoff,tp,tpr,fn,tn,tnr,fp,fpr)
    vec.list[[x]] <- vec
    x = x+1
  }
  
  # Create the data frame, append
  df <- as.data.frame(do.call(rbind, vec.list))
  colnames(df) <-  c("model","prSize","bgSize","cutoff","tp","tpr","fn","tn","tnr","fp","fpr")
  ls[[i]] <- df
  i = i+1
}

# Create the final data frame
accmeas <- bind_rows(ls)
glimpse(accmeas)

# Calculate additional metrics (precision, recall, and F1)
accmeas <- accmeas %>%
  mutate(
    model = factor(model),
    accuracy = (tp+tn)/(tp+fp+fn+tn),
    precision = tp / (tp+fp),
    recall = tp / (tp+fn),
    f1 = 2*(precision*recall/(precision+recall)),
    mcc = (tp*tn-fp*fn) / sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  )
glimpse(accmeas)

# save out
write_csv(accmeas,"data/tabular/mod/results/best_model/southern_rockies_accmeas.csv")

# Clean up
rm(abs,df,ls,pres,test,test.part,vec.list)
gc()


############################################################
# Identify the optimum cutoff threshold for classification #
############################################################

# Initialize a data frame to store the results
opt_thresh <- data.frame(
  model=integer(), 
  cutoff_f1=numeric(), 
  cutoff_mcc=numeric(),
  f1=numeric(), 
  mcc=numeric()
)

for(m in models) {
  
  print(paste0("model: ",m))
  df <- accmeas %>% filter(model == m)
  
  # Find the optimal threshold based on F1 or MCC
  opt_f1 <- df %>% arrange(desc(f1)) %>% slice(1)
  opt_mcc <- df %>% arrange(desc(mcc)) %>% slice(1)
  
  # Add to the data frame
  opt_thresh <- rbind(
    opt_thresh, 
    data.frame(
      model = m, 
      cutoff_f1 = opt_f1$cutoff, 
      cutoff_mcc = opt_mcc$cutoff,
      f1 = opt_f1$f1,
      mcc = opt_mcc$mcc)
    )
}

glimpse(opt_thresh)

# Plot the cutoff values by the two metrics

# Reshape the data
opt_thresh_l <- opt_thresh %>%
  gather(key = "metric", value = "cutoff", -model, -f1, -mcc) %>% 
  filter(metric %in% c("cutoff_f1", "cutoff_mcc"))

# Create the box plot
ggplot(opt_thresh_l, aes(x = metric, y = cutoff, fill = metric)) +
  geom_boxplot() +
  labs(title = "Distribution of Optimal Cutoffs by Metric",
       x = "Metric",
       y = "Cutoff Value") +
  scale_fill_discrete(name = "Metric", labels = c("F1 Score", "MCC")) +
  theme_minimal()

# Find the optimum cutoff as the mean across folds
print(paste("Average Optimal Threshold based on F1:", round(mean(opt_thresh$cutoff_f1),4)))
print(paste("Average Optimal Threshold based on MCC:", round(mean(opt_thresh$cutoff_mcc),4)))

# Write to a file
write_csv(opt_thresh,'data/tabular/mod/results/best_model/southern_rockies_opt_thresh.csv')


#########################################################
# Accuracy of the reference data (LFEVT, ITSP, TreeMap) #
#########################################################

# Function to calculate accuracy and F1 score
calculate_metrics <- function(df, true_label, pred_label) {
  tp <- sum((df[[true_label]] == 1) & (df[[pred_label]] == 1))
  tn <- sum((df[[true_label]] == 0) & (df[[pred_label]] == 0))
  fp <- sum((df[[true_label]] == 0) & (df[[pred_label]] == 1))
  fn <- sum((df[[true_label]] == 1) & (df[[pred_label]] == 0))
  
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  
  return(list(accuracy = accuracy, f1_score = f1_score, precision = precision, recall = recall))
}

# For each reference set
refs = c('lfevt','itsp','treemap')
# Get the optimal cutoff (mean based on F1 across folds)
opt_cutoff <- mean(opt_thresh$cutoff_f1)

# Reload the test data,
# Assign the classification label based on the optimum cutoff
testData <- read_csv(testPart) %>%
  rename(gee_id = `system:index`,
         TrueLabel = label) %>%
  dplyr::select(-.geo) %>%
  mutate(probability = probability*0.001, # scale back
         ClassLabel = if_else(probability >= opt_cutoff, 1, 0)) 
glimpse(testData)

# First, get the overall F1/MCC for the classification
overall_metrics <- calculate_metrics(testData, "TrueLabel", "ClassLabel")
print(overall_metrics)

# Initialize data frame to store results
ref_results_df <- data.frame(reference = character(), metric = character(), value = numeric())
# Apply the function to each reference column and to TrueLabel
for (ref in c(refs, 'TrueLabel')) {
  # Compare to predicted ClassLabel
  metrics_pred <- calculate_metrics(testData, ref, "ClassLabel")
  metrics_pred$reference <- ref
  metrics_pred$metric <- 'Predicted'
  
  # Compare to TrueLabel (actual label from training)
  metrics_true <- calculate_metrics(testData, ref, "TrueLabel")
  metrics_true$reference <- ref
  metrics_true$metric <- 'True'
  
  # Combine and add to the results data frame
  ref_results_df <- rbind(ref_results_df, metrics_pred, metrics_true)
}

# Convert reference and metric columns to factors
ref_results_df <- ref_results_df %>%
  mutate(reference = factor(reference),
         metric = factor(metric))

# Print and view the results data frame
print(ref_results_df)

# Write out to a file
write_csv(ref_results_df,'data/tabular/mod/results/best_model/southern_rockies_ref_accmeas.csv')

# Clean up
rm(list=ls())
gc()

