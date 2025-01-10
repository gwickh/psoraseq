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

class GPTrainer:
    def __init__(self, data_path, window_size=100):
        self.data_path = data_path
        self.window_size = window_size
        self.data = pd.read_csv(data_path)
        self.x_values = np.array(self.data.iloc[:, 0], dtype=int)
        self.y_values = self._preprocess_y_values(np.array(self.data.iloc[:, 1], dtype=np.float32))
        self.x_binned = np.floor(np.arange(len(self.y_values)) / window_size) * window_size

    def _preprocess_y_values(self, y):
        y = np.nan_to_num(y, nan=np.nan, posinf=np.nan, neginf=np.nan)
        y = pd.Series(y).replace(0, np.nan)
        y_interpolated = y.interpolate(method='spline', order=3).to_numpy()
        y_rolling = pd.Series(y_interpolated).rolling(
            window=self.window_size, 
            center=True,
            min_periods=1
        ).mean().to_numpy()
        return y_rolling

    def train_model(self):
        with pm.Model() as self.model:
            variance = pm.HalfNormal('variance', sigma=1)
            trend_kernel = pm.gp.cov.ExpQuad(1, ls=10)
            periodic_kernel = pm.gp.cov.Periodic(1, period=5.0, ls=10)
            cov = variance * (trend_kernel + periodic_kernel)
            gp = pm.gp.Latent(cov_func=cov)

            f = gp.prior("f", X=self.x_binned[:, None])

            sigma = pm.HalfNormal("sigma", sigma=100)
            nu = 1 + pm.Gamma("nu", alpha=2, beta=0.1)
            y_ = pm.StudentT("y", mu=f, lam=1.0 / sigma, nu=nu, observed=np.log2(self.y_values))

            self.approx = pm.fit(n=1000, method='advi')

    def plot_results(self):
        from pymc.gp.util import plot_gp_dist

        fig = plt.figure(figsize=(10, 4))
        ax = fig.gca()
        f_post = az.extract(self.approx, var_names="f").transpose("sample", ...)
        plot_gp_dist(ax, f_post, self.x_values[:, None])

        ax.plot(self.x_values[:, None], self.y_values, "ok", ms=3, label="Observed data")

        plt.xlabel("X")
        plt.ylabel("True f(x)")
        plt.title("Posterior distribution over $f(x)$ at the observed values")
        plt.legend()
        plt.show()

if __name__ == "__main__":
    trainer = GPTrainer(data_path='../../../relative_change.csv')
    trainer.train_model()
    trainer.plot_results()
