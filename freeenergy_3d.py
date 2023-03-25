import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from mpl_toolkits.mplot3d import Axes3D
plt.rcParams['font.size']=12

T = 310.15 #温度
kb = 1.38e-23 #ボルツマン定数

a = 20 #1辺
A = a*a #格子
Ns = 100 #actin分子数
n0 = 1e-3

#nf_values = [-10000*n0, -7500*n0, -5000*n0, -2500*n0, n0]

#color = ['green', 'red', 'blue', 'pink', 'black']

x_range = np.linspace(0.1, 0.9, 10)
nf_range = np.linspace(-100000*n0-n0, n0-n0, 10)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def S_func(x, A, Ns):
    return -((2*Ns-Ns*x)*np.log(Ns)/(2*A) + Ns*(1-x)*np.log(1-x)/A + Ns*x*np.log(x)/(2*A) + (A-Ns+Ns*x)*np.log(A-Ns+Ns*x)/A + (1-Ns*x/(2*A))*np.log(2*A-Ns*x) + (A-2*Ns+3*Ns*x)*np.log(x)/A + Ns*x*np.log(2)/A)


def E_func(x, n0, nf, Na):
    return n0 * Ns + Na * (nf - n0) * x

G = np.zeros((len(x_range), len(nf_range)))
for i, nf in enumerate(nf_range):
    for j, x in enumerate(x_range):
        E = E_func(x, n0, nf, Ns)
        S = S_func(x, A, Ns)
        G[j,i] = E - T * S

x, nf = np.meshgrid(x_range, nf_range)
ax.plot_surface(x, nf, G, cmap = "summer", alpha=0.5)
ax.plot(x.ravel(),nf.ravel(),G.ravel(),'co', markersize=3)

ax.set_title("G=E-TS")
ax.set_xlabel("x")
ax.set_ylabel("nf-n0")
ax.set_zlabel("G")
ax.legend()

plt.show()
