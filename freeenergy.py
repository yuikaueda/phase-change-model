import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

T = 310.15 #温度
kb = 1.38e-23 #ボルツマン定数

a = 20 #1辺
A = a*a #格子
Ns = 100 #actin分子数
n0 = 1e-3

nf_values = [-10000*n0, -7500*n0, -5000*n0, -2500*n0, n0]

color = ['green', 'red', 'blue', 'pink', 'black']

x_range = np.linspace(0.1, 0.9, 10000)

def S_func(x, A, Ns):
    return -((2*Ns-Ns*x)*np.log(Ns)/(2*A) + Ns*(1-x)*np.log(1-x)/A + Ns*x*np.log(x)/(2*A) + (A-Ns+Ns*x)*np.log(A-Ns+Ns*x)/A + (1-Ns*x/(2*A))*np.log(2*A-Ns*x) + (A-2*Ns+3*Ns*x)*np.log(x)/A + Ns*x*np.log(2)/A)

def E_func(x, n0, nf, Na):
    return n0*Ns + Na*(nf-n0)*x

for i in range(len(nf_values)):
    y = E_func(x_range, n0, nf_values[i], Ns) - T*S_func(x_range, A, Ns)
   # y =  -T*S_func(x_range, A, Ns)
    #y = E_func(x_range, n0, nf_values[i], Ns) 
    plt.plot(x_range, y, color=color[i], label=f"nf={nf_values[i]}")


plt.title("G=E-TS")
plt.xlabel("x")
plt.ylabel("G")
plt.legend()

# グラフを表示する
plt.show()
