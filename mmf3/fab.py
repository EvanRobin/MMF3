import numpy as np
import math

def exp_rek(x, epsilon): #b
    "Vrijednost eksponencijalne funkcije negativnog eksponenta dobivena rayvojem u red."
    result = 1.0
    e = epsilon
    k = 1
    mem = 1.0
    while abs(mem) > e:
        mem *= -x/k
        result += mem
        k += 1
    return result, k

def exp_red(x, epsilon): #a
    "Vrijednost eksponencijalne funkcije negativnog eksponenta dobivena rekurzivnom formulom."
    result = 1.0
    e = epsilon
    k = 1
    mem = 1.0
    while abs(mem) > e:
        mem = ((-1)**k)*(x**k)/math.factorial(k)
        result += mem
        k += 1    
    return result, k

def exp_rec(x, epsilon): #c
    "Vrijednost eksponencijalne funkcije negativnog eksponenta dobivena reciprociranjem standardne eksponencijalne funkcije iz a."
    result = 1.0
    e = epsilon
    k = 1
    mem = 1.0
    while abs(mem) > e:
        mem = (x**k)/math.factorial(k)
        result += mem
        k += 1
    return 1/result, k



limit = 1e-10

for i in range(0, 10):
    i *= 10
    print("%13.5e %13.5e %13.5e %d" %(exp_red(i, limit)[0], exp_rek(i, limit)[0], exp_rec(i, limit)[0], i))

print('\n%13.5e %13.5e %13.5e %d' %(exp_red(20, limit)[0], exp_rek(20, limit)[0], exp_rec(20, limit)[0], 20))