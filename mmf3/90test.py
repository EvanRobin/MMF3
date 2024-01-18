'''
 9. [3+5+2] Izraˇ cunajte vjerojatnost da atom neona mase m = 3.37 · 10−26 kg u idealnom plinu na temperaturi
 T = 300K ima brzinu manju ili jednaku srednjoj brzini vsrednje = 8kBT
 πm = 559.4m/s. Funkcija gusto´ ce
 vjerojatnosti za iznos brzina odgovara Maxwellovoj distribuciji
 3/2
 fMB(v) =
 m
 2πkBT
 2kBT
 4πv2e− mv2
 uz konstantu kB = 1.38064852 · 10−23 m2 kg s−2 K−1.
 a) Postavite problem i prilagodite ga za metodu numeriˇ cke integracije te napišite na papir.
 b) Gauss-Legendrovom metodom na¯ dite rješenja. Pokažite kod.
 c) Napišite rezultat na papir i broj toˇ caka integracije koji se koristio pri raˇ cunu. Komentirajte rezultat, koja je
 preciznost rezultata te kako ste je odredili.
'''
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

m = 3.37*10**(-26)
T = 300
v_av = 559.4
k = 1.38064852*10**(-23)

def probability_density_func(v):
    return ((m/(2*np.pi*k*T))**(3/2))*(4*np.pi*(v**2)*math.exp(-(m*v**2)/(2*k*T)))

print(int_gauleg(probability_density_func, 0, v_av, 500))
print(int_gauleg(probability_density_func, 0, v_av, 10))
