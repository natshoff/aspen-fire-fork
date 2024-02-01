"""
Plotting the results from Random Forest classification scenarios (different sets of input features)
For the best performing model, test for multicolinearity, do a grid search for hyperparameter tuning


"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score

import time
begin = time.time()

maindir = '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim1'

# Read back in the SBKCV results
results = pd.read_csv('data/tabular/mod/results/scenarios/rf_accmets_scenarios.csv')
# "feat_imps" includes feature importance values for each scenario
feat_imps = pd.read_csv('data/tabular/mod/results/scenarios/rf_feat_imps_scenarios.csv')
# "prob_preds" are the predicted probabilities for the reference data
prob_preds = pd.read_csv('data/tabular/mod/results/scenarios/rf_prob_preds_scenarios.csv')


###################################################################
# Metrics plots based on the default classification cutoff (0.50) #
###################################################################

# Plot the results across folds using boxplots

# Ensure 'Class' and 'Feature_Set' are categorical
results['Feature_Set'] = results['Feature_Set'].astype('category')

# Prep the accuracy metrics for facet wrapping
results_m = pd.melt(
    results,
    id_vars=['Fold', 'Feature_Set'],
    var_name='Metric',
    value_name='Score'
)

# Create a figure and axes
fig, axes = plt.subplots(1, 1, figsize=(6, 6))  # Adjust figure size as needed

# Plot MCC score across models
sns.boxplot(x='Feature_Set', y='Score', data=results_m[results_m['Metric'] == 'F1'], ax=axes)
axes.set_xlabel('Scenario')
axes.set_ylabel('F1')
axes.tick_params(axis='x', rotation=45)
# Adjust layout for readability
plt.tight_layout()
plt.subplots_adjust(top=0.9)  # Adjust top spacing for the title

plt.show()


############################################
# Select the best performing model (by F1) #
############################################

# Calculate mean across folds for each Feature_Set
mean_f1 = results.groupby('Feature_Set')['F1'].mean()  # or median()
# Display top models based on MCC and F1
print("Top models based on F1:\n", mean_f1.sort_values(ascending=False))

best_model_results = results[results['Feature_Set'] == mean_f1.idxmax()]
best_model_feat_imps = feat_imps[feat_imps['Feature_Set'] == mean_f1.idxmax()]
best_model_probs = prob_preds[prob_preds['Feature_Set'] == mean_f1.idxmax()]

# # Plot the best model feature importances
# plt.figure(figsize=(12, 8))  # Adjust the figure size as needed
# sns.boxplot(data=best_model_feat_imps, x='Feature', y='Importance')
#
# plt.title(f"Feature Importances for Best Model")
# plt.xticks(rotation=45, ha='right')  # Rotate feature names for better readability
# plt.ylabel('Importance')
# plt.xlabel('Feature')
#
# plt.tight_layout()
# plt.show()

# Write to file
best_model_results.to_csv('data/tabular/mod/results/scenarios/rf_results_best_model.csv')
best_model_feat_imps.to_csv('data/tabular/mod/results/scenarios/rf_feat_imps_best_model.csv')
best_model_probs.to_csv('data/tabular/mod/results/scenarios/rf_prob_preds_best_model.csv')


#####################
# ROC and PR curves #
#####################

# ROC Curve
fpr, tpr, roc_thresholds = roc_curve(best_model_probs['TrueLabel'], best_model_probs['PredictedProb'])
roc_auc = auc(fpr, tpr)
roc_data = pd.DataFrame({
    'False Positive Rate': fpr,
    'True Positive Rate': tpr,
    'Threshold': roc_thresholds
})
roc_data['ROC AUC'] = roc_auc  # Add a column for ROC AUC

# PR Curve
precision, recall, pr_thresholds = precision_recall_curve(
    best_model_probs['TrueLabel'],
    best_model_probs['PredictedProb']
)
pr_auc = average_precision_score(best_model_probs['TrueLabel'], best_model_probs['PredictedProb'])
pr_data = pd.DataFrame({
    'Precision': precision[:-1],  # Truncate the last value
    'Recall': recall[:-1],        # Truncate the last value
    'Threshold': pr_thresholds
})
pr_data['PR AUC'] = pr_auc  # Add a column for PR AUC

# Plot the ROC curve
plt.figure(figsize=(12, 6))

# First subplot: ROC Curve
plt.subplot(1, 2, 1)
plt.plot(roc_data['False Positive Rate'], roc_data['True Positive Rate'], color='blue', label=f'ROC curve (area = {roc_data["ROC AUC"].iloc[0]:.2f})')
plt.plot([0, 1], [0, 1], color='gray', linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve')
plt.legend(loc="lower right")

# Second subplot: PR Curve
plt.subplot(1, 2, 2)
plt.plot(pr_data['Recall'], pr_data['Precision'], color='green', label=f'PR curve (area = {pr_data["PR AUC"].iloc[0]:.2f})')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve')
plt.legend(loc="lower left")

plt.tight_layout()
plt.show()

# ####################################################################################################
# # Metrics plots based on moving cutoffs, incl. ROC & Prec-Rec curves for the best performing model #
# ####################################################################################################
#
# cutoffs = np.linspace(0, 1, num=10)  # 10 cutoffs from 0 to 1
# opt_cuts = pd.DataFrame()
# accmeas = pd.DataFrame()
#
# for ft_name, ft_data in ref_dict.items():
#     print(f'Running threshold optimization for reference set: {ft_name}')
#
#     sc_preds = prob_preds[(prob_preds['Feature_Set'] == ft_name)]
#
#     # Initialize a DataFrame to store the metrics for each cutoff per fold
#     fold_cutoffs = pd.DataFrame()
#
#     # Iterate through each fold
#     for fold in range(1, n_folds + 1):
#         # Filter the probabilities for the current class and fold
#         fold_probs = sc_preds[
#             (sc_preds['Fold'] == fold)
#         ]
#
#         # Initialize a list to store metrics for each cutoff for the current fold
#         fold_metrics = []
#
#         # Iterate through each cutoff
#         for c in cutoffs:
#             # Apply the cutoff
#             y_pred = (fold_probs['PredictedProb'] >= c).astype(int)
#
#             # Calculate metrics
#             f1 = f1_score(fold_probs['TrueLabel'], y_pred, zero_division=0)
#             mcc = matthews_corrcoef(fold_probs['TrueLabel'], y_pred)
#
#             # Append the metrics to the list
#             fold_metrics.append({
#                 'Feature_Set': ft_name,
#                 'Fold': fold,
#                 'Cutoff': c,
#                 'F1': f1,
#                 'MCC': mcc
#             })
#
#         # Convert the fold's metrics to a DataFrame and append to the fold_cutoff_metrics DataFrame
#         fold_metrics_df = pd.DataFrame(fold_metrics)
#         fold_cutoffs = pd.concat([fold_cutoffs, fold_metrics_df], ignore_index=True)
#
#     # Store the output for that scenario in a data frame
#     accmeas = pd.concat([accmeas,fold_cutoffs])
#
#     # Now fold_cutoffs contains all the metrics for each cutoff for each fold
#     # Group by cutoff and compute the mean of the metrics across folds
#     cutoff_avgs = fold_cutoffs.groupby('Cutoff').mean().reset_index()
#
#     # Find the cutoff with the highest average MCC (or F1)
#     opt_cuts = cutoff_avgs.loc[cutoff_avgs['MCC'].idxmax()]
#     opt_cuts['Feature_Set'] = ft_name
#
#     # Append the optimal cutoff for the current class to the optimal_cutoffs DataFrame
#     opt_cuts = pd.concat([
#         opt_cuts, pd.DataFrame([opt_cuts])], ignore_index=True
#     )
#
#
# print(opt_cuts)