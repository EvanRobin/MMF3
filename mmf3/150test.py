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
 15. [5+3+2] Sve su veliˇ cine dane u SI jedinicama. Koriste´ ci eksplicitnu shemu na¯ dite rješenja parcijalne
 diferencijalne jednadžbe
 ∂ρ(x,t)
 ∂t
 =0.02· ∂2ρ(x,t)
 ∂x2 .
 a) Rubne uvjete ρ(0,t) i ρ(20,t) uˇ citajte redom iz prvog i drugog stupca priložene datoteke (redni broj retka
 predstavlja trenutak j = t/0.5 > 0). Za poˇ cetni uvjet uzmite ρ(x,0) = −0.1x2 + 2x.
 b) Na istom grafu usporedite rješenja dobivena u trenucima t = j∆t za j = 150,300,450 i ∆t = 0.5.
 c) Koliki ∆x ima smisla uzeti? Što se doga¯ da s rubovima i šiljkom poˇ cetnog oblika gusto´ ce?
'''

D = 0.02

dt = 0.5
dx = 0.1
t = [0.0*dt, 150.0*dt, 300.0*dt, 400.0*dt]

def rho(x):
    if x >= 0.0 and x <= 20.0:
        return -0.1*x**2 + 2*x
    else:
        return 0.0



P1 = [0.0, 20.0, 0.0, t[0], dx, dt]
P2 = [0.0, 20.0, 0.0, t[1], dx, dt]
P3 = [0.0, 20.0, 0.0, t[2], dx, dt]
P4 = [0.0, 20.0, 0.0, t[3], dx, dt]
rub = [0.0, 0.0]

exp1 = dif(rho, P1, rub, D, 'exp')
exp2 = dif(rho, P2, rub, D, 'exp')
exp3 = dif(rho, P3, rub, D, 'exp')
exp4 = dif(rho, P4, rub, D, 'exp')

X = [x for x in np.arange(0.0, 20.0+dx, dx)]

fig = plt.figure(figsize=(7,5), dpi=120)
axes = fig.add_axes([0.10, 0.10, 0.85, 0.85])
plt.rcParams.update({'font.size': 8}) #type:ignore
axes.plot(X, exp1, lw=1.5, color='green')
axes.plot(X, exp2, lw=1.5, color='green', linestyle='-')
axes.plot(X, exp3, lw=1.5, color='green', linestyle='--')
axes.plot(X, exp4, lw=1.5, color='red')


axes.grid(lw=0.4)
axes.set_xlabel('x / m')
axes.set_ylabel('$\u03C1(x,t)$ / kgm$^{-1}$')
axes.legend(loc='best')
axes.set_title('Difuzija')
plt.show()

