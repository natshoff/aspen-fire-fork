"""
Tuning of aspen forest classification model for the Southern Rockies

    1. Load and prepare the presence (aspen) and background (non-aspen) learning data
    2. Create data frames for each of the classification scenarios
    3. Run an RF model for each classification scenario
        - Retrieve the feature importances,
        - Calculate the accuracy metrics (Prec/Rec, F1, MCC),
        - Test for class imbalance (?)
    4. Generate a boxplot of model performance for the different classification scenarios
    5. Identify the best performing model

"""

import os
import pandas as pd
import geopandas as gpd
import numpy as np

import time
begin = time.time()

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, matthews_corrcoef


maindir = '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim1'

# Load and prep the reference data

# Presence samples
Pres = os.path.join(maindir, 'data/spatial/mod/training/sampled/presence/pi_points_srme_m500_sampled.gpkg')
Pres = gpd.read_file(Pres)
Pres['PresAbs'] = 1
# Grab a list of the variable names
keep_cols = list(Pres.drop(['id','geometry'], axis=1).columns)
keep_cols.append('Block_ID')
X_cols = list(Pres.drop(['id','geometry','PresAbs'], axis=1).columns)

# Load the spatial blocks, with the count of presence data attached
blocks = gpd.read_file('data/spatial/mod/boundaries/spatial_block_grid_50km2_count.gpkg')
# Filter the presence points based on block statistics
# Need at least 100 presence points in the block
Pres = gpd.sjoin(Pres, blocks, how='left', predicate='within')
Pres = Pres.dropna(subset=['Block_ID'])
Pres = Pres[Pres.n_presence >= 100]
Pres = pd.DataFrame(Pres)

# Background samples
Bg = os.path.join(maindir, 'data/tabular/mod/training/background/bg_points_evt_sbcls_sampled.csv')
Bg = pd.read_csv(Bg)
Bg['PresAbs'] = 0
# Handle NULL values
Bg = Bg.dropna()

# Merge to create a presence/background learning dataset
pbl = pd.concat([Pres, Bg])
pbl = pbl[keep_cols]

# Check for NaNs in each column
nan_count = pbl.isna().sum()
print("NaN counts per column:\n", nan_count)

# Set the dependent and independent variables
y = pbl[['PresAbs']]  # 0 = non-aspen, 1 = aspen
# Retain numeric columns for model predictors
X = pbl[X_cols]

del keep_cols, blocks, Pres, Bg

# Check for the class imbalance ratio (should be ~10:1
len(pbl[pbl.PresAbs == 1]) / len(pbl[pbl.PresAbs == 0])


###############################################
# Set up the different classification scenarios

# Summer radar
sc1 = pbl[['PresAbs','Block_ID','VV_summer','VH_summer']]
sc2 = pbl[['PresAbs','Block_ID','VV_summer','VH_summer',
           'VH_summer_contrast','VH_summer_corr','VH_summer_ent','VH_summer_var',
           'VV_summer_contrast','VV_summer_corr', 'VV_summer_ent', 'VV_summer_var']]
# Winter radar
sc3 = pbl[['PresAbs','Block_ID','VV_winter','VH_winter']]
sc4 = pbl[['PresAbs','Block_ID','VV_winter','VH_winter',
           'VH_winter_contrast', 'VH_winter_corr','VH_winter_ent', 'VH_winter_var',
           'VV_winter_contrast', 'VV_winter_corr', 'VV_winter_ent','VV_winter_var']]
# Winter+Summer combined
sc5 = pbl[['PresAbs','Block_ID','VV_summer','VH_summer','VV_winter','VH_winter',
           'VH_summer_contrast', 'VH_summer_corr', 'VH_summer_ent', 'VH_summer_var',
           'VV_summer_contrast', 'VV_summer_corr', 'VV_summer_ent', 'VV_summer_var',
           'VH_winter_contrast', 'VH_winter_corr','VH_winter_ent', 'VH_winter_var',
           'VV_winter_contrast', 'VV_winter_corr', 'VV_winter_ent','VV_winter_var']]


# Summer spectral
sc6 = pbl[['PresAbs','Block_ID','B2_summer','B3_summer','B4_summer','B5_summer','B6_summer',
           'B7_summer','B8A_summer','B8_summer','B11_summer','B12_summer']]
sc7 = pbl[['PresAbs','Block_ID','B2_summer','B3_summer','B4_summer','B5_summer','B6_summer',
           'B7_summer','B8A_summer','B8_summer','B11_summer','B12_summer',
           'CIRE_summer','IRECI_summer','MCARI_summer','MNDWI_summer',
           'NDVI705_summer','SLAVI_summer']]

# Autumn spectral
sc8 = pbl[['PresAbs','Block_ID','B2_autumn','B3_autumn','B4_autumn','B5_autumn','B6_autumn',
           'B7_autumn','B8A_autumn','B8_autumn','B11_autumn','B12_autumn']]
sc9 = pbl[['PresAbs','Block_ID','B2_autumn','B3_autumn','B4_autumn','B5_autumn','B6_autumn',
           'B7_autumn','B8A_autumn','B8_autumn','B11_autumn','B12_autumn',
           'CIRE_autumn','IRECI_autumn','MCARI_autumn','MNDWI_autumn',
           'NDVI705_autumn','SLAVI_autumn']]

# Summer+Autumn Spectral
sc10 = pbl[['PresAbs','Block_ID','B2_summer','B3_summer','B4_summer','B5_summer','B6_summer',
            'B7_summer','B8A_summer','B8_summer','B11_summer','B12_summer',
            'CIRE_summer','IRECI_summer','MCARI_summer','MNDWI_summer','NDVI705_summer','SLAVI_summer',
            'B2_autumn', 'B3_autumn', 'B4_autumn', 'B5_autumn', 'B6_autumn',
            'B7_autumn', 'B8A_autumn', 'B8_autumn', 'B11_autumn', 'B12_autumn',
            'CIRE_autumn', 'IRECI_autumn', 'MCARI_autumn', 'MNDWI_autumn',
            'NDVI705_autumn', 'SLAVI_autumn']]


# Set up a dictionary to hold the classification scenarios
# Combine the reference datasets into a dictionary
ref_dict = {
    'Summer_S1': sc1,
    'Summer_S1_GLCM': sc2,
    'Winter_S1': sc3,
    'Winter_S1_GLCM': sc4,
    'Summer_Winter_S1': sc5,
    'Summer_S2': sc6,
    'Summer_S2_VI': sc7,
    'Autumn_S2': sc8,
    'Autumn_S2_VI': sc9,
    'Summer_Autumn_S2': sc10,
    'Combined_S1_S2': pbl  # This is just the full sampled data
}

del sc1,sc2,sc3,sc4,sc5,sc6,sc7,sc8,sc9,sc10


######################################
# Set up the classification workflow #

# Grab unique block IDs for spatial block cross validation
block_ids = pbl['Block_ID'].unique()

n_folds = 10  # Number of folds to run for SBKCV

# Define dataframes to store results for this feature set
results = pd.DataFrame()  # to store the model performance metrics
feat_imps = pd.DataFrame()  # to store the feature importances
prob_preds = pd.DataFrame()  # for testing optimum cutoff

# Loop the dictionary and run the classification models
for ft_name, ft_data in ref_dict.items():
    print(f'Running Random Forest for reference set: {ft_name}')

    # Loop through each fold
    for fold_idx in range(1, n_folds + 1):
        print(f'Running for fold {fold_idx} ...')

        # Split blocks into training and testing sets (70:30 split) with a different random state each time
        train_blocks, test_blocks = train_test_split(block_ids, test_size=0.3, random_state=fold_idx)

        # Set up the model data
        y_train = ft_data[ft_data['Block_ID'].isin(train_blocks)]['PresAbs']
        y_test = ft_data[ft_data['Block_ID'].isin(test_blocks)]['PresAbs']
        X_train = ft_data[ft_data['Block_ID'].isin(train_blocks)].drop(['PresAbs', 'Block_ID'], axis=1)
        X_test = ft_data[ft_data['Block_ID'].isin(test_blocks)].drop(['PresAbs', 'Block_ID'], axis=1)

        # Count presence and background samples in training and testing sets
        n_pres_train = np.sum(y_train == 1)
        n_bg_train = np.sum(y_train == 0)
        n_pres_test = np.sum(y_test == 1)
        n_bg_test = np.sum(y_test == 0)

        # For the 'scale_pos_weight' handling class imbalance
        wt = np.sum(y == 0) / np.sum(y == 1)  # currently unused

        # Initialize the Random Forest classifier
        rf_model = RandomForestClassifier(
            n_estimators=101,
            random_state=42,
            class_weight="balanced"
        )

        # Fit the model
        rf_model.fit(X_train, y_train)

        # Store feature importance
        fold_imps = pd.DataFrame({
            'Feature_Set': ft_name,
            'Fold': fold_idx,
            'Feature': X_train.columns,
            'Importance': rf_model.feature_importances_})

        feat_imps = pd.concat([feat_imps, fold_imps], axis=0)

        # Predict on the test set
        y_pred = rf_model.predict(X_test)

        # Retrieve the accuracy/performance metrics
        accuracy = accuracy_score(y_test, y_pred)
        precision = precision_score(y_test, y_pred, zero_division=0)
        recall = recall_score(y_test, y_pred, zero_division=0)
        f1 = f1_score(y_test, y_pred)
        mcc = matthews_corrcoef(y_test, y_pred)

        # Store the metrics into the results data frame
        fold_results = pd.DataFrame({
            'Feature_Set': [ft_name],
            'Fold': [fold_idx],
            'Accuracy': [accuracy],
            'Precision': [precision],
            'F1': [f1],
            'MCC': [mcc],
            'N_Pres_Train': [n_pres_train],
            'N_Bg_Train': [n_bg_train],
            'N_Pres_Test': [n_pres_test],
            'N_Bg_Test': [n_bg_test]
        })
        results = pd.concat([results, fold_results], ignore_index=True)

        # Store the probability values for cutoff testing
        y_pred_proba = rf_model.predict_proba(X_test)[:, 1]

        # Store probabilities and true labels
        fold_probs = pd.DataFrame({
            'Feature_Set': ft_name,
            'TrueLabel': y_test,
            'PredictedProb': y_pred_proba,
            'Fold': fold_idx
        })
        prob_preds = pd.concat([prob_preds, fold_probs], ignore_index=True)

        fold_idx += 1

    print(f"Time to complete classification scenario: {(time.time() - begin) / 60} minutes.")

print(f"Total elapsed time: {(time.time() - begin)/60} minutes.")

# Examine the results and export to tables

# "results" includes the accuracy assessment across folds for the different classification scenarios
print(results)
results.to_csv('data/tabular/mod/results/scenarios/rf_accmets_scenarios.csv')

# "feat_imps" includes feature importance values for each scenario
print(feat_imps)
feat_imps.to_csv('data/tabular/mod/results/scenarios/rf_feat_imps_scenarios.csv')

# "prob_preds" are the predicted probabilities for the reference data
print(prob_preds)
prob_preds.to_csv('data/tabular/mod/results/scenarios/rf_prob_preds_scenarios.csv')