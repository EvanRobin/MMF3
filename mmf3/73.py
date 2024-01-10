import math
import numpy as np
import matplotlib.pyplot as plt
import sys

def Euler(t0, x0, v0, function, N, tN):
    h = (tN-t0)/N 
    T = [t0]
    X = [x0]
    V = [v0]
    A = [function(t0, x0, v0)]
    while T[-1] <= tN:
        T.append(T[-1]+h)
        X.append(X[-1]+V[-1]*h)
        V.append(V[-1]+A[-1]*h)
        A.append(function(T[-1], X[-1], V[-1]))
    return T, X, V, A

def Runge_Kutta_2(t0, x0, v0, function, N, tN):
    h = (tN-t0)/N
    T = [t0]
    X = [x0]
    V = [v0]
    A = [function(t0, x0, v0)]
    while T[-1] <= tN:
        t = T[-1]
        x = X[-1]
        v = V[-1]
        k1v = function(t, x, v)
        k1x = v
        k2v = function(t+h, x+h*k1x, v+h*k1v)
        k2x = v+h*k1v
        T.append(T[-1]+h)
        X.append(x+(k1x+k2x)*h/2)
        V.append(v+(k1v+k2v)*h/2)
        A.append(function(T[-1], X[-1], V[-1]))
    return T, X, V, A

def Runge_Kutta_4(t0, x0, v0, function, N, tN):
    h = (tN-t0)/N
    T = [t0]
    X = [x0]
    V = [v0]
    A = [function(t0, x0, v0)]
    while T[-1] <= tN:
        t = T[-1]
        x = X[-1]
        v = V[-1]
        k1v = function(t, x, v)
        k1x = v
        k2v = function(t+h/2, x+h/2*k1x, v+h/2*k1v)
        k2x = v+h/2*k1v
        k3v = function(t+h/2, x+h/2*k2x, v+h/2*k2v)
        k3x = v+h/2*k2v
        k4v = function(t+h, x+h*k3x, v+h*k3v)
        k4x = v+h*k3v
        T.append(t+h)
        X.append(x+(k1x+2*k2x+2*k3x+k4x)*h/6)
        V.append(v+(k1v+2*k2v+2*k3v+k4v)*h/6)
        A.append(function(T[-1], X[-1], V[-1]))
    return T, X, V, A

def Euler_n(X0, function, N, xN):
    n = len(X0)
    h = (xN-X0[0])/N
    X = np.zeros((n+1, N+1))
    X[:-1, 0] = X0
    X[-1, 0] = function(X0)
    for i in range(N):
        X[0, i+1] = (i+1)*h
        X[1:-1, i+1] = X[1:-1, i]+h*X[2:, i]
        X[-1, i+1] = function(X[:-1, i])
    return X

def JUG(t0, X0, V0, acceleration, N, tN):
    h = (tN-t0)/N
    X = np.zeros((3, N+1))
    V = np.zeros((3, N+1))
    A = np.zeros((3, N+1))
    T = np.zeros(N+1)
    X[:, 0] = X0
    V[:, 0] = V0
    A[:, 0] = acceleration(t0, X0, V0)
    T[0] = t0
    for i in range(len(X)):
        for j in range(N):
            V[i, j+1] = V[i, j] + h*A[i, j]
            X[i, j+1] = X[i, j] + h*V[i, j] + 0.5*(h**2)*A[i, j]
            A[i, j+1] = acceleration(T[j], X[:, j], V[:, j])[i]
            T[j+1] = (j+1)*h
    return T, X, V, A

sys.getdefaultencoding()

X0_1D = [0.0, 4*np.pi/180, 0.0]
X0_3D = [4*np.pi/180, 0.0, 0.0]
V0_3D = [0.0, 0.0, 0.0]
l = 0.2484902028828339
m = 0.2
g = 9.81
N = 20000

T = 2*np.pi*np.sqrt(l/g)
t0 = 0.0
tN = 20*T

def fun_1D(vrijednosti):
    return -g/(l*m)*np.sin(vrijednosti[1])

def fun_Runge_Kutta_4(t, x, v):
    return -g/(l*m)*np.sin(x)

def fun_3D(t, x_vri, v_vri):
    return [-g/(l*m)*np.sin(x_vri[0]), 0.0, 0.0]

X_En , V_En = Euler_n(X0_1D, fun_1D, N, tN)[1], Euler_n(X0_1D, fun_1D, N, tN)[2]
X_RK, V_RK = Runge_Kutta_4(t0, X0_1D[1], X0_1D[2], fun_Runge_Kutta_4, N, tN)[1], Runge_Kutta_4(t0, X0_1D[1], X0_1D[2], fun_Runge_Kutta_4, N, tN)[2]
X_JUG, V_JUG = JUG(t0, X0_3D, V0_3D, fun_3D, N, tN)[1][0], JUG(t0, X0_3D, V0_3D, fun_3D, N, tN)[2][0]

fig = plt.figure(figsize=(7,6), dpi=120)
axes = fig.add_axes([0.15, 0.15, 0.75, 0.70])
plt.rcParams.update({'font.size': 8})
axes.plot(X_En, V_En, color='blue')
axes.plot(X_JUG, V_JUG, color='red')
axes.plot(X_RK, V_RK, color='green')
axes.grid(lw=0.5)
plt.show()