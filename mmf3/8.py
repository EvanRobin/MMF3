import numpy as np
import matplotlib.pyplot as plt 
import sys

def D_exp(g_x, xt0_N, D):
    dx = xt0_N[4]
    dt = xt0_N[5]
    N = int(xt0_N[1]/dx)
    M = int(xt0_N[3]/dt)
    alpha = D*dt/(dx**2)
    dif_p = np.zeros(N+1)
    dif_r = np.zeros(N+1)
    for u in range(len(dif_p)):
        dif_p[u] = g_x(xt0_N[0]+u*dx)
    for j in range(M+1):
        for i in range(1, N):
            dif_r[i] = alpha*dif_p[i+1]+(1-2*alpha)*dif_p[i]+alpha*dif_p[i-1]
        dif_p = np.copy(dif_r)
    return dif_r

sys.getdefaultencoding()

def rho(x):
    if x >= 2.0 and x <= 5.0:
        return 5.5
    else:
        return 0.0

dx, dt = 0.2, 0.5

j = [0, 100, 200, 300, 400]
t = [dt*J for J in j]

P1 = [0.0, 20.0, 0.0, t[0], dx, dt]
P2 = [0.0, 20.0, 0.0, t[1], dx, dt]
P3 = [0.0, 20.0, 0.0, t[2], dx, dt]
P4 = [0.0, 20.0, 0.0, t[3], dx, dt]
P5 = [0.0, 20.0, 0.0, t[4], dx, dt]

Dif = 1e-2

D1 = D_exp(rho, P1, Dif)
D2 = D_exp(rho, P2, Dif)
D3 = D_exp(rho, P3, Dif)
D4 = D_exp(rho, P4, Dif)
D5 = D_exp(rho, P5, Dif)

X = [x/dx for x in np.arange(0.0, 20.0+dx, dx)]

fig = plt.figure(figsize=(7,5), dpi=120)
axes = fig.add_axes([0.10, 0.10, 0.85, 0.85])
plt.rcParams.update({'font.size': 8})
axes.plot(X, D1, color='purple')
axes.plot(X, D2, color='blue')
axes.plot(X, D3, color='cyan')
axes.plot(X, D4, color='green')
axes.plot(X, D5, color='orange')
plt.show()