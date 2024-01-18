import math
import numpy as np
import matplotlib.pyplot as plt

def Runge_Kutta(t0, x0, v0, function, N, tN):
    h = (tN-t0)/N
    T = [t0]
    X = [x0]
    V = [v0]
    A = [function(t0, x0, v0)]
    F = [function(t0, x0, v0)*m]
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
        F.append(function(T[-1], X[-1], V[-1])*m)
    return T, X, V, A, F

'''
 10. [1+5+3+1] Na tijelo mase m = 2 kg, koje u poˇ cetnom trenutku miruje u ishodištu, djeluje ukupna sila
 F =2kgs−2·x−5kgs−1·v+0.5N.
 a) Napišite diferencijalne jednadžbe koje opisuju spomenuto gibanje i prilagodite ih za upotrebu RK4 metode.
 b) Simulirajte prvih 5 s opisano gibanje primjenom RK4 metode. Pokažite kod.
 c) Prikažite na istom grafu kako se položaj, brzina i sila mijenjaju u vremenu.
 d) Zapišite na papir konaˇ cne vrijednosti vremenskog koraka, položaja i brzine.
'''
x0= 0.0
v0 = 0.0
m = 2
N = 6000
tN = 5

def Force(x,v):
    return 2*x - 5*v + 0.5

def Accel_Startdust_synchron(t, x, v):
    return Force(x, v)/m


t, x, v, a, f = Runge_Kutta(0, x0, v0, Accel_Startdust_synchron, N, tN)[0], Runge_Kutta(0, x0, v0, Accel_Startdust_synchron, N, tN)[1], Runge_Kutta(0, x0, v0, Accel_Startdust_synchron, N, tN)[2], Runge_Kutta(0, x0, v0, Accel_Startdust_synchron, N, tN)[3], Runge_Kutta(0, x0, v0, Accel_Startdust_synchron, N, tN)[4]

fig = plt.figure(figsize=(11,5), dpi=110)
axes = fig.add_axes([0.15, 0.15, 0.75, 0.70])
plt.rcParams.update({'font.size': 9})
axes.plot(t, x, color='blue')
axes.plot(t, v, color='green')
axes.plot(t, a, color='red')
axes.plot(t, f, color='purple')
plt.show()
