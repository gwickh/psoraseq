import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pymc as pm


import numpy as np
import matplotlib.pyplot as plt

# Define parameters for the sinusoidal data
period = 100  # Single period from 0 to 100
num_points = 50  # 50 evenly spaced intervals
x = np.linspace(0, period, num_points)  # Generate x values
true_signal = np.sin(2 * np.pi * x / period)  # Sinusoidal function
noise = np.random.normal(0, 0.5, size=num_points)  # Generate noise with larger variance
y = true_signal + noise  # Add noise to the signal

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(x, true_signal, label="True Signal (Sinusoidal)", color="blue", linewidth=2)
plt.scatter(x, y, color="red", label="Noisy Data", alpha=0.7)
plt.title("Noisy Sinusoidal Data")
plt.xlabel("X")
plt.ylabel("Y")
plt.legend()
plt.grid(True)
plt.show()

####################################

with pm.Model() as model:
    variance = pm.HalfNormal('variance', sigma=1)
    trend_kernel = pm.gp.cov.ExpQuad(1, ls=10)  # Trend kernel
    periodic_kernel = pm.gp.cov.Periodic(1, period=5.0, ls=10)  # Periodic kernel
    cov = variance * (trend_kernel + periodic_kernel)
    gp = pm.gp.Latent(cov_func=cov)

    f = gp.prior("f", X=x[:, None])

    sigma = pm.HalfNormal("sigma", sigma=100)
    nu = 1 + pm.Gamma(
        "nu", alpha=2, beta=0.1
    )
    y_ = pm.StudentT("y", mu=f, lam=1.0 / sigma, nu=nu, observed=y)

    idata = pm.sample(1000, return_inferencedata=True, tune=500)
    
    
#####################################

# plot the samples from the gp posterior with samples and shading
from pymc.gp.util import plot_gp_dist

fig = plt.figure(figsize=(10, 4))
ax = fig.gca()
f_post = az.extract(idata, var_names="f").transpose("sample", ...)
plot_gp_dist(ax, f_post, x)

# plot the data and the true latent function
ax.plot(x, y, "ok", ms=3, label="Observed data")

# axis labels and title
plt.xlabel("X")
plt.ylabel("True f(x)")
plt.title("Posterior distribution over $f(x)$ at the observed values")
plt.legend()
