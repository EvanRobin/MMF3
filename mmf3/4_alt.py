import numpy as np
'''
a = np.array([[4, -2, 1],[3, 6, -4],[2, 1, 8]])
b = np.array([20, -30, 40])
x = np.linalg.solve(a, b)
'''
D = 7
a = [0, -2, -2, -2]#in the first eqation there is no a that is why its a 0
b = [4, 6, 6, 4]
c = [-2, -2, -2, -2]

d = [5, 0, 0, 0]

#xj+ej*xj+1=fj
#ej+1 = (cj+1)/(bj+1 - aj+1*ej)
#fj+1 = (bj+1 - aj+1*fj)/(bj+1 - aj+1*ej)
#all for (an*en-1 - bn)*xn = an*fn-1 - bn

e = [0, 0, 0, 0]
f = [0, 0, 0, 0]
x = [0, 0, 0, 0]

e[0] = c[0] / b[0]
f[0] = d[0] / b[0]

for i in range(1,D-1):#it will go from [1] to [D-1]
    den = (b[i] - a[i]*e[i-1]) #its the same in both formulas
    e[i] = c[i] / den
    f[i] = (d[i] - a[i]*f[i-1]) / den

x[D-1] = (a[D-1]*f[D-2] - d[D-1]) / (a[D-1]*e[D-2] - b[D-1])
for j in range(D-2, -1, -1):#goes till 0
    x[j] = f[j] - e[j]*x[j+1]
print('x=', [round(l,D+1) for l in x])