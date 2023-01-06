# Shortcut keyboard
## shift + enter: run selection/line in Python terminal

# Import packages
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import make_column_transformer
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.metrics import accuracy_score, confusion_matrix, precision_score, recall_score, roc_auc_score, roc_curve, f1_score
import itertools

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

### For LR5y
# Split the data into training and test sets
seed = 989 # Make the reporducible results
X_train, X_test, y_train, y_test = train_test_split(X, y_lr, test_size=0.3, random_state=seed)

# Make missing values as a separate category
X_train = X_train.fillna('missing')
X_test = X_test.fillna('missing')
y_train = y_train.fillna('missing')
y_test = y_test.fillna('missing')

# Recode the categorical features - new way
# Make a list of categorical variables to encode
X_train.dtypes
features_to_encode = X_train.columns[X_train.dtypes==object].tolist()

# Make a constructor to handle categorical features
col_trans = make_column_transformer(
    (OneHotEncoder(),features_to_encode),
    remainder = "passthrough")

# Train the random forest classifier
rf_classifier = RandomForestClassifier(
    min_samples_leaf=50,
    n_estimators=150,
    bootstrap=True,
    oob_score=True, # Use out-of-bag samples to estimate the generalization accuracy
    n_jobs=-1,
    random_state=seed,
    max_features='auto')

# Combine the constructor and classifier
rf_pipe = make_pipeline(col_trans, rf_classifier)
rf_pipe.fit(X_train, y_train) 
# rf_pipe includes 2 components: 
# 1. A constructor to handle inputs with categorical variables and transform into a correct type
# 2. A classifier that receives those newly transformed inputs from the constructor

# Make predictions on the test set
y_pred = rf_pipe.predict(X_test)

# Evaluate the classifer's performance
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

# Plot the ROC curve
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

# Confusion matrix
def plot_confusion_matrix(cm, classes, normalize = False,
                          title='Confusion matrix',
                          cmap=plt.cm.Greens): # can change color 
    plt.figure(figsize = (10, 10))
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title, size = 24)
    plt.colorbar(aspect=4)
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45, size = 14)
    plt.yticks(tick_marks, classes, size = 14)
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    # Label the plot
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt), fontsize = 20, horizontalalignment="center", color="white" if cm[i, j] > thresh else "black")
    plt.grid(None)
    plt.tight_layout()
    plt.ylabel('True label', size = 18)
    plt.xlabel('Predicted label', size = 18)

cm = confusion_matrix(y_test, y_pred)
plot_confusion_matrix(cm, classes = ['0 - No local recurrence', '1 - Yes local recurrence'], title = 'Local recurrence status within 5 years confusion matrix')

