import math
import numpy as np

def exp_series(x,epsilon):
    k=1
    result = 1.0
    member = 1.0
    while abs(member) > epsilon:
        member *= -x/k
        result += member
        k += 1
    return k,result

def exp_recursion(x, epsilon): 
    k = 1
    result = 1.0
    member = 1.0
    while abs(member) > epsilon:
        member = ((-1)**k)*(x**k)/math.factorial(k)
        result += member
        k += 1    
    return k, result

def exp_to_the_power_of(x, epsilon): 
    k = 1
    member = 1.0
    result = 1.0
    while abs(member) > epsilon:
        member = (x**k)/math.factorial(k)
        result += member
        k += 1
    return k, 1/result
eps=10**(-10)
for i in range (0,100,10):
    print(exp_series(i,eps),exp_recursion(i,eps),exp_to_the_power_of(i,eps))