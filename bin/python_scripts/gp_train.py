#!/usr/bin/env python3
# Filename:         gp_train.py
# Author:           Gregory Wickham
# Created:          2025-01-07
# Version:          1.2
# Date modified:    2025-01-08
# Description:      Train Latent Guassian Process on data

import pymc as pm
import arviz as az
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import CubicSpline

data = pd.read_csv('../../../relative_change.csv')

#get values
x_values = np.array(data.iloc[:, 0], dtype=int)[::4]
x_normalized = (x_values - min(x_values)) / (max(x_values) - min(x_values))

y = np.array(data.iloc[:, 1], dtype=np.float32)
y = np.nan_to_num(y, nan=np.nan, posinf=np.nan, neginf=np.nan)
y = pd.Series(y)
y = y.replace(0, np.nan)
y_interpolated = y.interpolate(method='spline', order=3).to_numpy()[::4]


#Create latent GP model    
with pm.Model() as model:
    variance = pm.HalfNormal('variance', sigma=1)
    trend_kernel = pm.gp.cov.ExpQuad(1, ls=10)  # Trend kernel
    periodic_kernel = pm.gp.cov.Periodic(1, period=5.0, ls=10)  # Periodic kernel
    cov = variance * (trend_kernel + periodic_kernel)
    gp = pm.gp.Latent(cov_func=cov)

    f = gp.prior("f", X=x_values[:, None])

    sigma = pm.HalfNormal("sigma", sigma=100)
    nu = 1 + pm.Gamma(
        "nu", alpha=2, beta=0.1
    )
    y_ = pm.StudentT("y", mu=f, lam=1.0 / sigma, nu=nu, observed=y_interpolated)

    trace = pm.sample(1000, chains=4, tune=500, target_accept=0.9, max_treedepth=10)