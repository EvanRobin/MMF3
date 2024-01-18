import math
import numpy as np
import matplotlib.pyplot as plt

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

''' 
14. [2+6+1+1] Koriste´ ci RK4 metodu riješite diferencijalnu jednadžbu koja opisuje titranje tijela mase m =
 20 g pod djelovanjem sile ⃗ F = −k⃗x − b⃗v. Za konstantu opruge uzmite k = 80 mN m −1, a za konstantu
 gušenja b = 4 g s−1. Neka tijelo u poˇ cetnom trenutku miruje na udaljenosti 5 cm od ishodišta.
 a) Napišite diferencijalnu jednadžbu za položaj danog tijela i njeno analitiˇ cko rješenje.
 b) Usporedite numeriˇ cko i analitiˇ cko rješenje tijekom prvih 20 s tako da ispišete po stupcima redom vrijeme,
 položaj odre¯ den analitiˇ cki i numeriˇ cki. Prepišite usporedbe za zadnji lokalni maksimum i odabrani ∆t.
 c) Što se doga¯ da s amplitudom oscilacija tijekom vremena?
 d) Odredite ∆t za koji je u zadnjem uspore¯ denom lokalnom maksimumu relativna greška manja od 0.5%?
''' 

x0= 0.05
v0 = 0.0
m = 0.02
k = 80
b = 0.004
N = 6000
tN = 20
w = math.sqrt(k/m)

def Force(x,v):
    return - k*x - b*v

def Accel_Startdust_synchron(t, x, v):
    return Force(x, v)/m

def ANAL_osci(t, x):
    return 0.05*math.exp(-(b/(2*m))*t)*np.cos(w*t + 0)

def ANAL(t0, x0, v0,function, N, tN):
    h = (tN-t0)/N
    T = [t0]
    X = [x0]
    while T[-1] <= tN:
        t = T[-1]
        x = X[-1]
        T.append(t+h)
        X.append(function(T[-1], X[-1]))
    return T, X

who_cares = 0

t_ANAL, x_ANAL = ANAL(0, x0, v0, ANAL_osci, N, tN)[0], ANAL(0, x0, v0, ANAL_osci, N, tN)[1]

t, x = Runge_Kutta(0, x0, v0, Accel_Startdust_synchron, N, tN)[0], Runge_Kutta(0, x0, v0, Accel_Startdust_synchron, N, tN)[1]
fig = plt.figure(figsize=(11,5), dpi=110)
axes = fig.add_axes([0.15, 0.15, 0.75, 0.70])
plt.rcParams.update({'font.size': 9})
axes.plot(t, x, color='blue')
axes.plot(t_ANAL, x_ANAL, color='red')

plt.show()