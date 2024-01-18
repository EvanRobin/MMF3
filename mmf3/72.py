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

def Runge_Kutta2(t0, x0, v0, function, N, tN):
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

def Euler_n(X0, function, N, xN):
    '''
    Eulerova metoda za rjesavanje obicnih diferencijalnih jednadzbi n-tog
    reda za jednodimenzionalni problem.
    \n
    \nX0 ---------- vektor pocetnih uvjeta
    \nfunction ---- funkcija derivacije n-tog reda
    \nN ----------- broj intervala
    \nxN ---------- krajnja vrijednost
    '''
    n = len(X0) #red diferencijalne jednadzbe
    h = (xN-X0[0])/N #korak
    X = np.zeros((n+1, N+1)) #matrica rjesenja
    X[:-1, 0] = X0 #dodavanje pocetnih uvjeta
    X[-1, 0] = function(X0) #dodavanje vrijednosti derivacije n-tog reda
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

X0 = [0.0, 0.0, 0.0]
m = 4.7
F = 5.2
k = 400
def Full_Metal_Alchamist(Y):
    return F/m
N = 6000
tN = 3.2

t, x, v, a = Euler_n(X0, Full_Metal_Alchamist, N, tN)[0], Euler_n(X0, Full_Metal_Alchamist, N, tN)[1], Euler_n(X0, Full_Metal_Alchamist, N, tN)[2], Euler_n(X0, Full_Metal_Alchamist, N, tN)[3]

fig = plt.figure(figsize=(11,5), dpi=110)
axes = fig.add_axes([0.15, 0.15, 0.75, 0.70])
plt.rcParams.update({'font.size': 9})
axes.plot(t, x, color='blue')
axes.plot(t, v, color='green')
axes.plot(t, a, color='red')
plt.show()