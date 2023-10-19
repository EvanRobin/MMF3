import numpy as np

a = np.array([[4, -2, 1],[3, 6, -4],[2, 1, 8]])
b = np.array([20, -30, 40])
x = np.linalg.solve(a, b)

print(x)
