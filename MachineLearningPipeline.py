# Import packages
import pandas as pd 
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

# Load the data
df = pd.read_csv('../data/derived/2022Dec21_dat_TNBC.csv')

# Manipulate the data
df.select = df[['Age.at.Diagnosis', 
                'LR5y', 'DR5y', 'NP5y']]
df.select = df.select.dropna()
dat = pd.DataFrame(df.select)

### Local recurrence
# Split the data into features and target
X.LR = dat.drop(['LR5y', 'DR5y', 'NP5y'], axis=1)
y.LR = dat['LR5y']

# Split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X.LR, y.LR, test_size=0.2, random_state=999)

# Train the model
model = RandomForestClassifier()
model.fit(X_train, y_train)

# Make predictions on the test set
predictions = model.predict(X_test)

# Evaluate the model's performance
accuracy = model.score(X_test, y_test)
print('Accuracy:', accuracy)

### New primary breast cancer
# Split the data into features and target
X.NP = dat.drop(['LR5y', 'DR5y', 'NP5y'], axis=1)
y.NP = dat['NP5y']

# Split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X.NP, y.NP, test_size=0.2, random_state=999)

# Train the model
model = RandomForestClassifier()
model.fit(X_train, y_train)

# Make predictions on the test set
predictions = model.predict(X_test)

# Evaluate the model's performance
accuracy = model.score(X_test, y_test)
print('Accuracy:', accuracy)
