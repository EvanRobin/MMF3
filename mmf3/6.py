import numpy as np

def integral_trapezoidal(f, a, b, m):
    h = (b - a)/m  #h se odabire
    f_a = f(a)  
    f_b = f(b)  #krajnji vrijednosti
    rez = (f_a + f_b)/2
    k = 1 #int koji koji dok je m manji zbraja trapezoide
    while k < m:
        rez += f(a + k*h)  #trap se zbroji u sumi njih
        k += 1
    return rez*h #pomnoži se sa h i dobije se pov

def integral_Simpson(f, a, b, m):
    h = (b - a)/m  
    f_a = f(a)  
    f_b = f(b)  
    rez = f_a + f_b
    k = 1
    while k < m:
        if k % 2 == 0:
            rez += 2*f(a + k*h)  #pol ako je poz parabola
        else:
            rez += 4*f(a + k*h) #pol ako je neg parabola
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

o = [10, 50, 100, 1000, 2000]

for i in range(len(o)):
    print((o[i], integral_trapezoidal(p, v - dv, v + dv, o[i]), integral_Simpson(p, v - dv, v + dv, o[i])))

h = 0.1
v = [1e3, 5e3, 1e4, 5e4]

for i in range(len(v)):
    print((v[i]/h, v[i], integral_trapezoidal(p, 0., v[i], v[i]/h), integral_Simpson(p, 0., v[i], v[i]/h)))