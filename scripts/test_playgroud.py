import warnings; warnings.filterwarnings("ignore")
import numpy as np
import matplotlib.pyplot as plt
import pcsaft
from pcsaft.plotting import figure, finish, save

for name in ["methane", "decane"]:
    c = pcsaft.get(name)
    d = pcsaft.segment_diameter(300.0, c)
    print(f"{name:8s} d/sigma = {d/c.sigma:.4f}")

# now plot the whole curve against reduced temperature
for name in ["methane", "decane"]:
    c = pcsaft.get(name)
    T = np.linspace(20, 2000, 300)
    d = np.array([pcsaft.segment_diameter(t, c) for t in T])
    plt.plot(T / c.eps_k, d / c.sigma, label=name)
plt.xlabel("T / (eps/k)"); plt.ylabel("d / sigma"); plt.legend(); plt.show()


fig, ax = figure()                      # sized for a LaTeX column
ax.plot(T / c.eps_k, d / c.sigma, label="something")        # ax is an ordinary matplotlib axis
finish(ax, "T / (eps/k)", "d / sigma", legend=True)
print(save(fig, "ex01_diameter"))       # -> figures/ex01_diameter.pdf and .png