import argparse
import os
import sys
import numpy as np
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.linear_model import ElasticNet
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor
import pandas as pd
import math
from math import sqrt
import sklearn.preprocessing as sk
from sklearn import metrics
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split
import random
from random import randint
from sklearn.model_selection import StratifiedKFold
from scipy.stats import pearsonr, spearmanr, kendalltau
from sklearn.metrics import mean_squared_error
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def RF(X_train, Y_train, X_test, n_estimators, depth, mtry, cv, seed):
    Y_pred_0 = np.zeros([X_test.shape[0], 1])
    RF_pipe= Pipeline([('scaler', StandardScaler()),('RnF', RandomForestRegressor())])
    model = GridSearchCV(RF_pipe, 
                         param_grid={"RnF__n_estimators": n_estimators, "RnF__max_depth": depth,
                                     "RnF__max_features": mtry},
                         scoring='neg_mean_squared_error', cv=KFold(n_splits=cv, shuffle=True, random_state=seed))
    y_train = Y_train
    x_train = X_train
    model.fit(x_train, y_train)
    Y_pred_0 = model.predict(X_test)
    return Y_pred_0


def EN(X_train, Y_train, X_test, alpha, l1ratio, cv, seed):
    Y_pred_0 = np.zeros([X_test.shape[0], 1])
    EN_pipe = Pipeline([('scaler', StandardScaler()),('En', ElasticNet())])
    model = GridSearchCV(EN_pipe, param_grid={"En__alpha": alpha, "En__l1_ratio": l1ratio}, 
                         scoring='neg_mean_squared_error', cv=KFold(n_splits=cv, shuffle=True, random_state=seed))
    y_train = Y_train
    x_train = X_train
    model.fit(x_train, y_train)
    Y_pred_0 = model.predict(X_test)
    return Y_pred_0

def RG(X_train, Y_train, X_test, alpha, cv, seed):
    Y_pred_0 = np.zeros([X_test.shape[0], 1])
    RG_pipe = Pipeline([('scaler', StandardScaler()),('RdG', Ridge())])
    model = GridSearchCV(RG_pipe, param_grid={"RdG__alpha": alpha}, 
                         scoring='neg_mean_squared_error', cv=KFold(n_splits=cv, shuffle=True, random_state=seed))
    y_train = Y_train
    x_train = X_train
    model.fit(x_train, y_train)
    Y_pred_0 = model.predict(X_test)
    return Y_pred_0