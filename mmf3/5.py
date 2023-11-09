import matplotlib.pyplot as plt
import numpy as np
import math

x1 = [1, 2, 3]
y1 = [-5, -12, -15]
x2 = [1, 2, 3, 4, 5]
y2 = [-5, -12, -15, -8, 15]

def l(x, y, t):
    p = 0
    for i in range(0,len(x)):
        l=1
        for j in range(0,len(x)):
            if j != i:
                l *= (t-x[j])/(x[i]-x[j])
        p += l*y[i]
    return p

t = np.arange(0.5, 5.5, 0.001)
p1 = []
p2 = []
for i in range (0,len(t)):
    p2.append(l(x2, y2, t[i]))
    p1.append(l(x1, y1, t[i]))

plt.plot(t, p2, c='green')
plt.plot(t, p1, c='blue') 
plt.scatter(x2, y2, c='red') 
plt.xlabel('x - axis')  
plt.ylabel('y - axis')  
plt.title('<(")')  
plt.show()  

#thats the first part
