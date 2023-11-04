
# LIBRARIES

library(tidyverse)

# DATA

sdir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim1/data/tabular/mod/results/sensitivity'

dfs.test <- list.files(sdir, pattern="^srme_skcv_test.*?\\.csv")

for(df in dfs.test){
  sdir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim1/data/tabular/mod/results/sensitivity'
  
  # Get the path
  path <- paste(sdir,df,sep = "/")
  print(path)
  
  # Get name
  name <- str_sub(df, start = 22, end = -5)
  print(name)
  
  # Read the test data
  test <- read_csv(path) %>%
    rename(gee_id = `system:index`) %>%
    dplyr::select(-.geo) %>%
    mutate(probability = probability*0.001) # scale back

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
    df$subset <- name
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
  write_csv(accmeas,paste0(sdir,'/accmeas/',name,'.csv'))
  rm(list = ls())
  gc()
}


# Now create a single data frame with the model averages for each subset of features

sdir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim1/data/tabular/mod/results/sensitivity'
adir <- paste0(sdir,'/accmeas/')
dfs.accmeas <- list.files(adir, pattern="*.csv")

i=1
ls <- list()
for(df in dfs.accmeas){
  
  # Read in the data frame
  path <- paste0(adir,df)
  print(path)
  
  accmeas <- read_csv(path)
  
  # Get name
  name <- str_sub(df, end = -5)
  print(name)
  
  # Calculate optimum threshold based on F1 statistic by model
  
  bestModel = accmeas %>%
    group_by(model) %>%
    summarize(f1Max = max(f1,na.rm=T),
              prMean = mean(prSize),
              bgMean = mean(bgSize))
  
  # Add the subset feature name as a column
  bestModel$subset <- name
  
  # Export to a list for binding
  ls[[i]] <- bestModel
  i = i+1
  
}

# Bind together the means
accmeas.means <- bind_rows(ls)
glimpse(accmeas.means)

# Write to a CSV
write_csv(accmeas.means,paste0(sdir,'/accmeas_f1best_sensitivity.csv'))
rm(list = ls())
gc()



