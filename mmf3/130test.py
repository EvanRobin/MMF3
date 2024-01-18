import numpy as np
import math

def gauleg(x_1, x_2, n):
    '''Metoda Gauss-Legendrove kvadrature za odredivanje nultocki x i njihovih statistickih
    tezina w. Granice integracije su x1 i x2, a duljina lista x i w je n.'''
    x = [0]*n #nul tocke 
    w = [0]*n #tezina nul tocke
    epsilon = 1e-6#min za nul tocku
    m = int((n+1)/2)#brojac
    x_m = (x_1+x_2)/2 #sredina interval 
    x_l = (x_2-x_1)/2 #pola interval
    for i in range(1, m+1):
        z = np.cos(np.pi*(i-0.25)/(n+0.5)) #poc nul tocka
        while True:
            #pocetni polinom sa L1 i L2
            L1 = 1.0
            L2 = 0.0
            for j in range(1, n+1):
                L3 = L2
                L2 = L1
                L1 = ((2*j-1)*z*L2-(j-1)*L3)/j
            dL = n*(z*L1-L2)/(z**2-1) #derivira se polinom
            z_0 = z
            z = z_0-L1/dL
            if abs(z-z_0) < epsilon:
                break
        x[i-1] = x_m-x_l*z #racuna se nul tocka
        x[n-i] = x_m+x_l*z
        w[i-1] = 2*x_l/((1-z**2)*dL**2) #izračunata težina
        w[n-i] = w[i-1] #simetrično je
    return x, w

def int_gauleg(function, a, b, m):
    x, w = gauleg(a, b, m)[0], gauleg(a, b, m)[1] #nul tocke i težine
    rez = 0.0
    for i in range(len(x)):
        rez += function(x[i])*w[i] #suma
    return rez 

def int_trapezoidal(function, a, b, m):
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

'''
 13. [10] Koriste´ ci trapeznu i Simpsonovu metodu odredite masu ravnog štapa ˇ ciji su poˇ cetak x = −0.5 i kraj
 x =2telinearna gusto´ ca λ = exp(x) + x5. Usporedite rezultate obiju metoda za isti broj koraka.
'''
x1 = -0.5
x2 = 2

def linear_distribution_of_mass(x):
    return math.exp(x) + x**5

print("Simp", int_Simpson(linear_distribution_of_mass, -0.5, 5, 100))
print("Trap", int_trapezoidal(linear_distribution_of_mass, -0.5, 5, 100))
print("Gauss-Legendrove", int_gauleg(linear_distribution_of_mass, -0.5, 5, 100))