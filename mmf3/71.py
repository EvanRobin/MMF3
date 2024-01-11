import math
import numpy as np
import matplotlib.pyplot as plt

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

def Prediktor_korektor(t0, x0, v0, function, N, tN):
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

def Runge_Kutta(t0, x0, v0, function, N, tN):
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

m = 0.2
l = 0.2484902028828339
y0 = 4
g = 9.81
v0 = 0.0

y0_r = y0*np.pi/180

def analiticko(t, x, v):
    return -g*x/l

def numericko(t, x, v):
    return -g/l*np.sin(x)

T = 2*np.pi*np.sqrt(l/g)
t1 = 1*T
t2 = 20*T


t = []
theta = []
for i in np.arange(t1, t2+1e-5, 1e-5):
    t.append(i)
    theta.append(y0_r*np.cos(i*np.sqrt(g/l)))

'''
fig = plt.figure(figsize=(11,5), dpi=110)
axes = fig.add_axes([0.10, 0.10, 0.80, 0.80])
plt.rcParams.update({'font.size': 9})           
axes.plot(Euler(0.0 , y0_r, v0, numericko, 100000, t2)[0], Euler(0.0 , y0_r, v0, numericko, 100000, t2)[1], color='blue', lw=1.5)
axes.plot(Euler(0.0 , y0_r, v0, numericko, 10000, t2)[0], Euler(0.0 , y0_r, v0, numericko, 10000, t2)[1], color='green', lw=1.5)
axes.plot(Runge_Kutta(0.0 , y0_r, v0, numericko, 100000, t2)[0], Runge_Kutta(0.0 , y0_r, v0, numericko, 100000, t2)[1], color='yellow')
axes.plot(Runge_Kutta(0.0 , y0_r, v0, numericko, 10000, t2)[0], Runge_Kutta(0.0 , y0_r, v0, numericko, 10000, t2)[1], color='red')
plt.show()
'''
'''
fig = plt.figure(figsize=(11,5), dpi=110)
axes = fig.add_axes([0.15, 0.15, 0.75, 0.70])
plt.rcParams.update({'font.size': 9})          
axes.plot(Euler(0.0 , y0_r, v0, numericko, 10000, t2)[0], Euler(0.0 , y0_r, v0, numericko, 10000, t2)[1], color='blue')
axes.plot(Euler(0.0 , y0_r, v0, numericko, 50000, t2)[0], Euler(0.0 , y0_r, v0, numericko, 50000, t2)[1], color='green')
axes.plot(Euler(0.0 , y0_r, v0, numericko, 100000, t2)[0], Euler(0.0 , y0_r, v0, numericko, 100000, t2)[1], color='orange')
axes.plot(t, theta, color='purple')
axes.set_xlim(t1, t2)
plt.show()
'''

''''''
t = []
theta_4 = []
theta_8 = []
theta_16 = []
theta_32 = []
for i in np.arange(0.0, 7*T+1e-5, 1e-5):
    t.append(i)
    theta_4.append(y0_r*np.cos(i*np.sqrt(g/l)))
    theta_8.append(y0_r*2*np.cos(i*np.sqrt(g/l)))
    theta_16.append(y0_r*4*np.cos(i*np.sqrt(g/l)))
    theta_32.append(y0_r*8*np.cos(i*np.sqrt(g/l)))

fig = plt.figure(figsize=(11,5), dpi=110)
axes = fig.add_axes([0.15, 0.15, 0.75, 0.70])
plt.rcParams.update({'font.size': 9})          
axes.plot(Runge_Kutta(0.0 , y0_r, v0, numericko, 1e5, 7*T)[0], Runge_Kutta(0.0 , y0_r, v0, numericko, 1e5, 7*T)[1], color='green')
axes.plot(Runge_Kutta(0.0 , y0_r*2, v0, numericko, 1e5, 7*T)[0], Runge_Kutta(0.0 , y0_r*2, v0, numericko, 1e5, 7*T)[1], color='red')
axes.plot(Runge_Kutta(0.0 , y0_r*4, v0, numericko, 1e5, 7*T)[0], Runge_Kutta(0.0 , y0_r*4, v0, numericko, 1e5, 7*T)[1], color='yellow')
axes.plot(Runge_Kutta(0.0 , y0_r*8, v0, numericko, 1e5, 7*T)[0], Runge_Kutta(0.0 , y0_r*8, v0, numericko, 1e5, 7*T)[1], color='blue', linestyle='--')
axes.plot(t, theta_32, color='blue')
plt.show()
