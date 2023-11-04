
# Libraries

library(tidyverse)

# Data

accmeas <- read_csv("data/tabular/mod/results/accmeas.csv")
glimpse(accmeas)

# Plot the ROC curve based on TPR and FPR

roc <- ggplot(data=accmeas) +
  geom_line(aes(x=fpr,y=tpr)) +
  geom_point(aes(x=fpr,y=tpr),size=0.3) +
  labs(x='False Positive Rate',y='True Positive Rate',title='ROC Curve') +
  theme(legend.position = "none")
roc

# Calculate optimum threshold base don GMean statistic

cutoffOpt = accmeas[which.max(accmeas$gmean),]$cutoff
fprOpt = accmeas[which.max(accmeas$gmean),]$fpr
tprOpt = accmeas[which.max(accmeas$gmean),]$tpr
  
# Add to our ROC plot

roc +
  geom_point(aes(x=fprOpt,y=tprOpt),color='#981220',size=4) +
  geom_text(aes(x=fprOpt,y=tprOpt),
            label = paste0('Optimal threshold \n for class: ', round(cutoffOpt,3)),
            nudge_x = 0.16, nudge_y = -0.12, size = 5) +
  labs(title='ROC Curve w/ Optimal Threshold Based on GMean Statistic') +
  theme_minimal(14)

# Create a Precision-Recall Curve

prc <- ggplot(data=accmeas) +
  geom_point(aes(x=recall,y=precision),size=0.4) +
  geom_line(aes(x=recall,y=precision)) +
  labs(x='Recall',y='Precision',title='Precision-Recall Curve') +
  theme_minimal(14)
prc

# Add in the optimal threshold based on the F1 score

cutoffOptF1 = accmeas[which.max(accmeas$f1),]$cutoff
precisionOptF1 = accmeas[which.max(accmeas$f1),]$precision
recallOptF1 = accmeas[which.max(accmeas$f1),]$recall

prc +
  geom_point(aes(x=recallOptF1,y=precisionOptF1),color='#981220',size=4) +
  geom_text(aes(x=recallOptF1,y=precisionOptF1),
            label = paste0('Optimal threshold \n for class: ', round(cutoffOptF1,3)),
            nudge_x = -0.20, nudge_y = -0.12, size = 5) +
  labs(title='ROC Curve w/ F1 Opt. Threshold') +
  theme_minimal(14)


# Calculate the accuracy metrics based on the two thresholds (GMean and F1-based)
# add to our plots

# GMean

print("GMean-based accuracy metrics:")
(f1Opt = accmeas[which.max(accmeas$gmean),]$f1)
(oaOpt = accmeas[which.max(accmeas$gmean),]$accuracy)
(precOpt = accmeas[which.max(accmeas$gmean),]$precision)
(recOpt = accmeas[which.max(accmeas$gmean),]$recall)

# F1

print("F1-based accuracy metrics:")
(f1Opt = accmeas[which.max(accmeas$f1),]$f1)
(oaOpt = accmeas[which.max(accmeas$f1),]$accuracy)
(precOpt = accmeas[which.max(accmeas$f1),]$precision)
(recOpt = accmeas[which.max(accmeas$f1),]$recall)


# Now, work with the average across model runs

means <- accmeas %>%
  mutate(cutoff = as.character(cutoff)) %>%
  group_by(cutoff) %>%
  summarise(
    fpr = mean(fpr,na.rm=T),
    tpr = mean(tpr,na.rm=T),
    precision = mean(precision,na.rm=T),
    recall = mean(recall,na.rm=T),
    gmean = mean(gmean,na.rm=T),
    f1 = mean(f1,na.rm=T),
    accuracy = mean(accuracy)
  ) %>%
  ungroup() %>%
  mutate(cutoff = round(as.double(cutoff),4))

# Recreate the plots

# Plot the ROC curve based on average TPR and FPR

roc <- ggplot(data=means) +
  geom_line(aes(x=fpr,y=tpr)) +
  geom_point(aes(x=fpr,y=tpr),size=0.6) +
  labs(x='False Positive Rate',y='True Positive Rate',title='ROC Curve') +
  theme(legend.position = "none")
roc

# Calculate optimum threshold base don GMean statistic

cutoffOpt = means[which.max(means$gmean),]$cutoff
fprOpt = means[which.max(means$gmean),]$fpr
tprOpt = means[which.max(means$gmean),]$tpr

# Add to our ROC plot

roc +
  geom_point(aes(x=fprOpt,y=tprOpt),color='#981220',size=4) +
  geom_text(aes(x=fprOpt,y=tprOpt),
            label = paste0('Optimal threshold \n for class: ', round(cutoffOpt,3)),
            nudge_x = 0.23, nudge_y = -0.12, size = 5) +
  labs(title='ROC Curve w/ GMean Opt. Threshold (means)') +
  theme_minimal(14)


# Create a Precision-Recall Curve

prc <- ggplot(data=means) +
  geom_point(aes(x=recall,y=precision),size=0.4) +
  geom_line(aes(x=recall,y=precision)) +
  labs(x='Recall',y='Precision',title='Precision-Recall Curve') +
  theme_minimal(14)
prc

# Add in the optimal threshold based on the F1 score

cutoffOptF1 = means[which.max(means$f1),]$cutoff
precisionOptF1 = means[which.max(means$f1),]$precision
recallOptF1 = means[which.max(means$f1),]$recall

prc +
  geom_point(aes(x=recallOptF1,y=precisionOptF1),color='#981220',size=4) +
  geom_text(aes(x=recallOptF1,y=precisionOptF1),
            label = paste0('Optimal threshold \n for class: ', round(cutoffOptF1,3)),
            nudge_x = -0.20, nudge_y = -0.12, size = 5) +
  labs(title='ROC Curve w/ F1 Opt. Threshold') +
  theme_minimal(14)

# add to our plots

# GMean

print("GMean-based accuracy metrics:")
(f1Opt = means[which.max(means$gmean),]$f1)
(oaOpt = means[which.max(means$gmean),]$accuracy)
(precOpt = means[which.max(means$gmean),]$precision)
(recOpt = means[which.max(means$gmean),]$recall)

# F1

print("F1-based accuracy metrics:")
(f1Opt = means[which.max(means$f1),]$f1)
(oaOpt = means[which.max(means$f1),]$accuracy)
(precOpt = means[which.max(means$f1),]$precision)
(recOpt = means[which.max(means$f1),]$recall)
