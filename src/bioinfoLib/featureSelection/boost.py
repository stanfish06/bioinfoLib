# Copyright 2025 Zhiyuan Yu (Heemskerk's lab, University of Michigan)

import numpy as np
from scipy.linalg import svd
from xgboost import XGBClassifier, XGBRegressor


def rank_features(X, y, mode="r"):
    if mode == "r":
        model = XGBRegressor()
    elif mode == "c":
        model = XGBClassifier()
    else:
        return None
    model.fit(X, y)
    return model.feature_importances_


if __name__ == "__main__":
    # test if XGBoost can interpret PCs or not
    np.random.seed(2)
    X = np.random.rand(100, 5)
    X = (X - np.mean(X, axis=0)) / np.std(X, axis=0)
    u, s, vt = svd(X)
    fs = rank_features(X, u[:, 0])
    print(fs / np.max(fs))
    print(np.abs(vt[0, :]) / np.max(np.abs(vt[0, :])))
