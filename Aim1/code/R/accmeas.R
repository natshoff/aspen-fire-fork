
# LIBRARIES

library(tidyverse)

# DATA

train <- "data/tabular/mod/results/prop/srme_skcv_train_probs_prop.csv"
test <- "data/tabular/mod/results/prop/srme_skcv_test_probs_prop.csv"

# Sampled training (probabilities)
train <- read_csv(train) %>%
  rename(gee_id = `system:index`) %>%
  dplyr::select(-.geo) %>%
  mutate(probability = probability*0.001) # scale back
glimpse(train)

# Sampled testing (probabilities)
test <- read_csv(test) %>%
  rename(gee_id = `system:index`) %>%
  dplyr::select(-.geo) %>%
  mutate(probability = probability*0.001) # scale back
glimpse(test)


# Calculate the accuracy metrics

models <- unique(test$seed)

# create a list to hold the results
i=1
ls <- list()
for(model in models){
  
  # subset to the seed #
  test.part <- test %>% filter(seed == model)
  
  # grab presence/background
  pres <- test.part %>% filter(label==1)
  abs <- test.part %>% filter(label==0)
  
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

# Calculate accuracy metrics
accmeas <- accmeas %>%
  mutate(
    model = factor(model),
    specificity = tn / (fp+tn),
    precision = tp / (tp+fp),
    recall = tp / (tp+fn),
    f1 = 2*(precision*recall/(precision+recall)),
    accuracy = (tp+tn)/(tp+fp+fn+tn),
    gmean = sqrt(tpr * (1-fpr)),
    mcc = (tp*tn-fp*fn) / sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  )

# save out
write_csv(accmeas,"data/tabular/mod/results/accmeas_prop.csv")
rm(list = ls())
gc()


#########################################################
# Accuracy of the reference data (LFEVT, ITSP, TreeMap) #
#########################################################

print("Calculating accuracy from the reference datasets ...")

ref.train <- read_csv('data/tabular/mod/results/reference/srme_skcv_mtry7/srme_ref_train_full_model_mtry7.csv')
ref.test <- read_csv('data/tabular/mod/results/reference/srme_skcv_mtry7/srme_ref_test_full_model_mtry7.csv')

glimpse(ref.test)

# Train data
ref.train.sum <- ref.train %>%
  group_by(seed) %>%
  summarize(
    lfevt_acc = sum(lfevt)/sum(label) * 100,
    itsp_acc = sum(itsp)/sum(label) * 100,
    treemap_acc = sum(treemap)/sum(label) * 100)
glimpse(ref.train.sum)
mean(ref.train.sum$lfevt_acc)
mean(ref.train.sum$itsp_acc)
mean(ref.train.sum$treemap_acc)

# Test data
ref.test.sum <- ref.test %>%
  group_by(seed) %>%
  summarize(
    lfevt_acc = sum(lfevt)/sum(label) * 100,
    itsp_acc = sum(itsp)/sum(label) * 100,
    treemap_acc = sum(treemap)/sum(label) * 100)
glimpse(ref.test.sum)
mean(ref.test.sum$lfevt_acc)
mean(ref.test.sum$itsp_acc)
mean(ref.test.sum$treemap_acc)


#############################################
# Accuracy of the independent validation data

val <- read_csv()

