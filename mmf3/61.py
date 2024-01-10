import numpy as np


def gauleg(x_1, x_2, n):
    '''Metoda Gauss-Legendrove kvadrature za odredivanje nultocki x i njihovih statistickih
    tezina w. Granice integracije su x1 i x2, a duljina lista x i w je n.'''
    x = [0]*n
    w = [0]*n
    epsilon = 1e-6
    m = int((n+1)/2)
    x_m = (x_1+x_2)/2
    x_l = (x_2-x_1)/2
    for i in range(1, m+1):
        z = np.cos(np.pi*(i-0.25)/(n+0.5)) 
        while True:
            L1 = 1.0
            L2 = 0.0
            for j in range(1, n+1):
                L3 = L2
                L2 = L1
                L1 = ((2*j-1)*z*L2-(j-1)*L3)/j
            dL = n*(z*L1-L2)/(z**2-1) 
            z_0 = z
            z = z_0-L1/dL
            if abs(z-z_0) < epsilon:
                break
        x[i-1] = x_m-x_l*z 
        x[n-i] = x_m+x_l*z
        w[i-1] = 2*x_l/((1-z**2)*dL**2) 
        w[n-i] = w[i-1] 
    return x, w

def int_gauleg(function, a, b, m):
    '''Metoda rjesavanja odredenog integrala koristeci Gauss-Legendrovu kvadraturu.'''
    x, w = gauleg(a, b, m)[0], gauleg(a, b, m)[1]
    rez = 0.0
    for i in range(len(x)):
        rez += function(x[i])*w[i]
    return rez 

def int_trapezoidal(function, a, b, m):
    '''Trapezna formula'''
    h = (b - a)/m 
    f_a = function(a)
    f_b = function(b)
    rez = (f_a + f_b)/2
    k = 1
    while k < m:
        rez += function(a + k*h)
        k += 1
    return rez*h


def int_Simpson(function, a, b, m):
    '''Simpsonova metoda'''
    h = (b - a)/m
    f_a = function(a)
    f_b = function(b)
    rez = f_a + f_b
    k = 1
    while k < m:
        if k % 2 == 0:
            rez += 2*function(a + k*h)
        else:
            rez += 4*function(a + k*h)
        k += 1
    return rez*h/3
  

m = 3.37e-26
T = 300
dv = 50
k = 1.38064852e-23

v = np.sqrt((8*k*T)/(np.pi*m))

K = m/(2*k*T)

def p(x):
    return (K/np.pi)**(3/2)*4*np.pi*(x**2)*np.exp(-K*(x**2))

o = [10, 50, 100]



for i in range(len(o)):
    print(o[i], int_trapezoidal(p, v-dv, v+dv, o[i]), int_Simpson(p, v - dv, v + dv, o[i]), int_gauleg(p, v-dv, v+dv, o[i])),