# Shortcut keyboard
## command + enter: run selection/line in Python terminal

# Import basic packages
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt

# Load the data
tnbc_df = pd.read_csv('../data/derived/2022Dec21_dat_TNBC.csv')
tnbc_df.shape # Dimension of data: (604, 115)

# Check missing values for all variables including targets
tnbc_df.isnull().sum()/len(tnbc_df)
tnbc_df.LR5y.isnull().sum()/len(tnbc_df)
tnbc_df.NP5y.isnull().sum()/len(tnbc_df)

# View all variable names
print(tnbc_df.columns.tolist())

# View the distribution of target and features
tnbc_df.LR5y.value_counts(dropna=False)
tnbc_df.NP5y.value_counts(dropna=False)
tnbc_df[['Clinical.T.Stage2']].value_counts(dropna=False)
tnbc_df[['Clinical.N.Stage2']].value_counts(dropna=False)

# Impute the missing values of target using the most frequent value
from sklearn.impute import SimpleImputer
simple_imp = SimpleImputer(missing_values = np.nan, strategy = "most_frequent")
simple_imp.fit(tnbc_df[['LR5y', 'NP5y']])
tnbc_df_filled = tnbc_df.copy()
tnbc_df_filled[['LR5y', 'NP5y']] = simple_imp.transform(tnbc_df_filled[['LR5y', 'NP5y']])
tnbc_df_filled.LR5y.value_counts(dropna=False)
tnbc_df_filled.NP5y.value_counts(dropna=False)

# Recode the categorical features - old way
# clinical_t_stage_map = {'T1|T2': 0, 'T3|T4': 1}
# clinical_n_stage_map = {'N0': 0, 'N1|N2|N3': 1}
# tnbc_df['Clinical.T.Stage2'] = tnbc_df['Clinical.T.Stage2'].map(clinical_t_stage_map)
# tnbc_df['Clinical.N.Stage2'] = tnbc_df['Clinical.N.Stage2'].map(clinical_n_stage_map)

# Manipulate the data
tnbc_df_select = tnbc_df_filled[['Age.at.Diagnosis', 'Clinical.T.Stage2', 'Clinical.N.Stage2', # features
                                 'LR5y', 'NP5y']] # target
dat = pd.DataFrame(tnbc_df_select)

# Split the data into features and target
X = dat.drop(['LR5y', 'NP5y'], axis=1)
y_lr = dat['LR5y'] 
y_np = dat['NP5y'] 

### For LR5y (5-year local recurrence)
# Split the data into training and test sets
from sklearn.model_selection import train_test_split
seed = 989 # Make the reporducible results
X_train, X_test, y_train, y_test = train_test_split(X, y_lr, test_size=0.3, random_state=seed)

# Make missing values as a separate category
X_train = X_train.fillna('missing')
X_test = X_test.fillna('missing')
y_train = y_train.fillna('missing')
y_test = y_test.fillna('missing')
# X_train[["Clinical.T.Stage2"]].value_counts(dropna=False)

# Recode the categorical features - new way
# Make a list of categorical variables to encode
X_train.dtypes
features_to_encode = X_train.columns[X_train.dtypes==object].tolist()

# Make a constructor to handle categorical features
from sklearn.compose import make_column_transformer
from sklearn.preprocessing import OneHotEncoder
col_trans = make_column_transformer(
    (OneHotEncoder(),features_to_encode),
    remainder = "passthrough")

# Train the random forest classifier
from sklearn.ensemble import RandomForestClassifier
rf_classifier = RandomForestClassifier(
    min_samples_leaf=50,
    n_estimators=150,
    bootstrap=True,
    oob_score=True, # Use out-of-bag samples to estimate the generalization accuracy
    n_jobs=-1,
    random_state=seed,
    max_features='auto')

# Combine the constructor and classifier
from sklearn.pipeline import make_pipeline
rf_pipe = make_pipeline(col_trans, rf_classifier)
rf_pipe.fit(X_train, y_train) 
# rf_pipe includes 2 components: 
# 1. A constructor to handle inputs with categorical variables and transform into a correct type
# 2. A classifier that receives those newly transformed inputs from the constructor

# Make predictions on the test set
y_pred = rf_pipe.predict(X_test)

# Evaluate the classifer's performance
from sklearn.metrics import accuracy_score, confusion_matrix, precision_score, recall_score, roc_auc_score, roc_curve, f1_score

# Accuracy = (TP + TN) / (TP+TN+FP+FN)
# Recall = TP / (TP + FN)
# Precision = TP / (TP + FP)
# F1-score = 2 * Precision * Recall / (Precision+Recall)

print(f'The accuracy of the classifer is {round(accuracy_score(y_test,y_pred)*100, 1)} %')

train_prob = rf_pipe.predict_proba(X_train)[:,1] 
# predict_proba(dataframe)[:,1] gives the predicted probability distribution of class label 1 from the dataframe
# The predicted probability is the fraction of samples of the same class label in a leaf
test_prob = rf_pipe.predict_proba(X_test)[:, 1]
train_pred = rf_pipe.predict(X_train)

print(f'Train ROC AUC Score: {roc_auc_score(y_train, train_prob)}')
print(f'Test ROC AUC  Score: {roc_auc_score(y_test, test_prob)}')

## Plot the ROC curve
def model_roc_curve(y_pred, test_prob, train_pred, train_prob):
    baseline = {}
    baseline['recall']=recall_score(y_test, [1 for _ in range(len(y_test))])
    baseline['precision'] = precision_score(y_test, [1 for _ in range(len(y_test))])
    baseline['roc'] = 0.5
    results = {}
    results['recall'] = recall_score(y_test, y_pred)
    results['precision'] = precision_score(y_test, y_pred)
    results['roc'] = roc_auc_score(y_test, test_prob)
    train_results = {}
    train_results['recall'] = recall_score(y_train, train_pred)
    train_results['precision'] = precision_score(y_train, train_pred)
    train_results['roc'] = roc_auc_score(y_train, train_prob)
    for metric in ['recall', 'precision', 'roc']:  
        print(f'{metric.capitalize()}, Baseline: {round(baseline[metric], 2)}, Test: {round(results[metric], 2)}, Train: {round(train_results[metric], 2)}')
     # Calculate false positive rates and true positive rates
    base_fpr, base_tpr, _ = roc_curve(y_test, [1 for _ in range(len(y_test))])
    model_fpr, model_tpr, _ = roc_curve(y_test, test_prob)
    plt.figure(figsize = (8, 6))
    plt.rcParams['font.size'] = 16
    # Plot both curves
    plt.plot(base_fpr, base_tpr, 'b', label = 'baseline')
    plt.plot(model_fpr, model_tpr, 'r', label = 'model')
    plt.legend();
    plt.xlabel('False Positive Rate');
    plt.ylabel('True Positive Rate'); plt.title('ROC Curves');
    plt.show();

model_roc_curve(y_pred, test_prob, train_pred, train_prob)

## Plot the confusion matrix
import itertools
def plot_confusion_matrix(cm, classes, normalize=False,
                          title = 'Confusion matrix',
                          cmap = plt.cm.Greens): # can change color 
    plt.figure(figsize = (10, 10))
    plt.imshow(cm, interpolation = 'nearest', cmap = cmap)
    plt.title(title, size = 24)
    plt.colorbar(aspect=4)
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45, size = 14)
    plt.yticks(tick_marks, classes, size = 14)
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    # Label the plot
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt), fontsize = 20, horizontalalignment = "center", color = "white" if cm[i, j] > thresh else "black")
        plt.grid(None)
        plt.tight_layout()
        plt.ylabel('True label', size = 18)
        plt.xlabel('Predicted label', size = 18)

cm = confusion_matrix(y_test, y_pred)
plot_confusion_matrix(cm, classes = ['0 - No LR', '1 - Yes LR'], title = '5-year local recurrence status confusion matrix')
plt.show()

## Estimate the feature importance
print(rf_classifier.feature_importances_)
print(f" There are {len(rf_classifier.feature_importances_)} features in total")

### Map values in the feature importance to the features in X_train
print(col_trans.fit_transform(X_train)[0,:]) # First rows of all columns
X_train.iloc[0,:] 
# The rf_classifier.feature_importances_ places all the encoded categorical variables before the numeric variables

### Create a proper encoded X_train
#### The function encode_and_bind encodes the categorical variables and then combine them with the original dataframe
def encode_and_bind(original_dataframe, features_to_encode):
    dummies = pd.get_dummies(original_dataframe[features_to_encode])
    res = pd.concat([dummies, original_dataframe], axis=1)
    res = res.drop(features_to_encode, axis=1)
    return(res)
X_train_encoded = encode_and_bind(X_train, features_to_encode)
X_train_encoded # [422 rows x 7 columns]

feature_importances = list(zip(X_train_encoded, rf_classifier.feature_importances_))
feature_importances
# [('Clinical.T.Stage2_T1|T2', 0.12097216088490814), ('Clinical.T.Stage2_T3|T4', 0.0), ('Clinical.T.Stage2_missing', 0.14133985245682396), ('Clinical.N.Stage2_N0', 0.2619941129812915), ('Clinical.N.Stage2_N1|N2|N3', 0.05964351032704749), ('Clinical.N.Stage2_missing', 0.10321231201103147), ('Age.at.Diagnosis', 0.31283805133889747)]

### Sort the feature importances by most important first
feature_importances_ranked = sorted(feature_importances, key = lambda x: x[1], reverse = True)

### Print out the feature and importances
[print('Feature: {:35} Importance: {}'.format(*pair)) for pair in feature_importances_ranked]

### Plot the top feature importance
feature_names_top = [i[0] for i in feature_importances_ranked[:3]]
y_ticks = np.arange(0, len(feature_names_top))
x_axis = [i[1] for i in feature_importances_ranked[:3]]
plt.figure(figsize = (6, 10))
plt.barh(feature_names_top, x_axis) # Make a horizontal barplot
plt.title('Random Forest Feature Importance (Top 3)',
          fontdict= {'fontname':'Comic Sans MS','fontsize' : 10})
plt.xlabel('Features', fontdict= {'fontsize' : 8})
plt.show()

# Tune the hyperparameters with RandomSearchCV
from pprint import pprint
print('Parameters currently in use:\n')
pprint(rf_classifier.get_params())
# {'bootstrap': True,
#  'ccp_alpha': 0.0,
#  'class_weight': None,
#  'criterion': 'gini',
#  'max_depth': None,
#  'max_features': 'auto',
#  'max_leaf_nodes': None,
#  'max_samples': None,
#  'min_impurity_decrease': 0.0,
#  'min_samples_leaf': 50,
#  'min_samples_split': 2,
#  'min_weight_fraction_leaf': 0.0,
#  'n_estimators': 150,
#  'n_jobs': -1,
#  'oob_score': True,
#  'random_state': 989,
#  'verbose': 0,
#  'warm_start': False}

## Create a grid of parameters for the model to randomly pick and train, hence the name Random Search.
from sklearn.model_selection import RandomizedSearchCV
n_estimators = [int(x) for x in np.linspace(start = 100, stop = 700, num = 50)]
# [100, 112, 124, 136, 148, 161, 173, 185, 197, 210, 222, 234, 246, 259, 271, 283, 295, 308, 320, 332, 344, 357, 369, 381, 393, 406, 418, 430, 442, 455, 467, 479, 491, 504, 516, 528, 540, 553, 565, 577, 589, 602, 614, 626, 638, 651, 663, 675, 687, 700]

