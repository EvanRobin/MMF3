import numpy as np
import matplotlib.pyplot as plt 
import sys

def Thomas(lower_diag, main_diag, upper_diag, solutions):
    '''Thomasova metoda za pronalazenje rjesenja linearnog sustava jednadzbi
    iz trodijagonalne matrice.'''
    a = lower_diag
    b = main_diag
    c = upper_diag
    d = solutions
    n = len(d)
    c_trenutno = [0]*n
    d_trenutno = [0]*n
    rez = [0]*n
    c_trenutno[0] = c[0]/b[0]
    d_trenutno[0] = d[0]/b[0]
    for i in range(1, n-1):
        c_trenutno[i] = c[i]/(b[i]-a[i-1]*c_trenutno[i-1])
        d_trenutno[i] = (d[i]-a[i-1]*d_trenutno[i-1])/(b[i]-a[i-1]*c_trenutno[i-1])
    d_trenutno[n-1] = (d[n-1]-a[n-2]*d_trenutno[n-2])/(b[n-1]-a[n-2]*c_trenutno[n-2])
    rez[n-1] = d_trenutno[n-1]
    for j in range(n-2, -1, -1):
        rez[j] = d_trenutno[j]-c_trenutno[j]*rez[j+1]
    return rez

def dif(g_x, xt0_N, R, D, met):
    '''Eksplicitna i implicitna metoda rjesavanje difuzijske parcijalne
    diferencijalne jednadzbe u 1D.
    \ng_x ------ funkcija pocetnih uvjeta za x varijablu
    \nxt0_N ---- vektor pocetnih uvjeta - [x0, xN, t0, tN, dx, dt]
    \nR -------- vektor rubnih uvjeta - [x<, x>]
    \nD -------- konstanta difuzije
    \nmet ------ metoda rjesavanja - "exp" / "imp"'''
    dx = xt0_N[4] #korak polozaja
    dt = xt0_N[5] #korak vremena
    N = int((xt0_N[1]-xt0_N[0])/dx) #broj tocaka u prostoru
    M = int((xt0_N[3]-xt0_N[2])/dt) #broj tocaka u vremenu
    dL = R[0] #lijevi rubni uvjet
    dD = R[1] #desni rubni uvjet
    alpha = D*dt/(dx**2)
    dif_p = np.zeros(N+1)
    dif_r = np.zeros(N+1)
    for u in range(len(dif_p)):
        dif_p[u] = g_x(xt0_N[0]+u*dx)
    if met == 'exp': #eksplicitna metoda
        for j in range(M+1): #vrijeme
            for i in range(1, N): #polozaj
                dif_r[i] = alpha*dif_p[i+1]+(1-2*alpha)*dif_p[i]+alpha*dif_p[i-1]
            dif_r[0], dif_r[-1] = dL, dD
            dif_p = np.copy(dif_r)
    elif met == 'imp': #implicitna metoda
        down = [-alpha]*N
        mid = [1+2*alpha]*(N+1)
        up = [-alpha]*N
        for j in range(M+1): #vrijeme
            dif_r = Thomas(down, mid, up, dif_p)
            dif_r[0], dif_r[-1] = dL, dD
            dif_p = np.copy(dif_r)
    else:
        print('Invalid method input.')
    return dif_r

sys.getdefaultencoding()

'''
11. [(3+3)+2+2] Linearna gusto´ ca tvari unutar 1D prostora opisana je jednadžbom
 ∂ρ(x,t)
 ∂t
 =1m2s−1 ∂2ρ(x,t)
 ∂x2 .
 (2.1)
 2
a) Prilagodite programe koji uz pomo´ c eksplicitne i implicitne sheme rješavaju (2.1) za sluˇ caja kada je
 ρ(±5 m,t) = 0 te ρ(x,0) = 25kgm−3·x2−1kgm−5·x4 ; x∈[−5m,5m]
 0.0 kg / m
 ;
 x / ∈ [−5 m,5 m] 
''' 
D = 1.0
t = [0.0, 4.0, 8.0]
dx = 0.1
dt = (dx**2)/2

def rho(x):
    if x >= -5.0 and x <= 5.0:
        return 25*(x*2)-(x*4)
    else:
        return 0.0

P0 = [-15.0, 15.0, 0.0, 0, dx, dt]
P4 = [-15.0, 15.0, 0.0, 4, dx, dt]
P8 = [-15.0, 15.0, 0.0, 8, dx, dt]
rub = [0.0, 0.0]

exp0 = dif(rho, P0, rub, D, 'exp')
exp4 = dif(rho, P4, rub, D, 'exp')
exp8 = dif(rho, P8, rub, D, 'exp')
imp0 = dif(rho, P0, rub, D, 'imp')
imp4 = dif(rho, P4, rub, D, 'imp')
imp8 = dif(rho, P8, rub, D, 'imp')

X = [x/dx for x in np.arange(-15.0, 15.0+dx, dx)]

'''
fig = plt.figure(figsize=(7,5), dpi=120)
axes = fig.add_axes([0.10, 0.10, 0.85, 0.85])
plt.rcParams.update({'font.size': 8}) #type:ignore
axes.plot(X, D1, label='t = {}$\u0394$t'.format(0.4), lw=1.4, color='purple')
axes.plot(X, D2, label='t = {}$\u0394$t'.format(8), lw=1.4, color='blue')
axes.grid(lw=0.5)
axes.set_xlabel('x / $\u0394$x')
axes.set_ylabel('$\u03C1(x,t)$ / kgm$^{-1}$')
axes.legend(loc='best')
axes.set_title('Difuzija')
plt.show()
'''
print(len(exp4),len(imp8),len(X))