import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.interpolate import CubicSpline

x = np.arange(10)
y = np.sin(x)

cs = CubicSpline(x, y)

xs = np.arange( -0, 10.1, 0.1)

fig, ax = plt.subplots(figsize = (6.5, 4))
ax.plot(x, y, 'o', label='data')
ax.plot(xs, np.sin(xs), label='True')
ax.plot(xs, cs(xs), label='S')
ax.plot(xs, cs(xs, 1), label="S'")
ax.plot(xs, cs(xs, 2), label="S''")
ax.plot(xs, cs(xs, 3), label="S'''")
ax.set_xlim(-0.5, 9.5)
ax.legend(loc='lower left', ncol=2)
plt.show()

print(x)