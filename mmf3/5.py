import matplotlib.pyplot as plt
import numpy as np
import math

x = [1, 2, 3, 4, 5]
y = [-5, -12, -15, -8, 15]

def l(x, y, t):
    p = 0
    for i in range(0,len(x)):
        l=1
        for j in range(0,len(x)):
            if j != i:
                l *= (t-x[j])/(x[i]-x[j])
        p += l*y[i]
    return p

t = np.arange(1, 5, 0.001)
p = np.array([])
for i in range (0,len(t)):
    append



