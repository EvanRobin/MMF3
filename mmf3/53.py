import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.interpolate import CubicSpline

f = open("v(H-H)_AK.txt", "a")
f.truncate(0)
cols = ("#", "r/A", "V/K")
f.write('\t'.join(cols) + '\n')
cols1 = ("#", "--------", "--------",)
f.write('\t'.join(cols1) + '\n')

def read_collum(file_name, convert=float, sep=None):
    with open(file_name) as file:
        for line in file:
            if line.startswith('#'):
                continue
            else:
                a = (convert(line.split(sep=sep)[0]))
                b = (convert(line.split(sep=sep)[1]))
            cols2 = (str(format((a*0.52917721092),'.6f')), str(format((b*315775.04), '.6f')))
            f.write('\t'.join(cols2) + '\n')

read_collum('V(H-H).txt')
f.close()

def write_list(file_name, convert=float, sep=None):
    Al = []
    Bl = []
    with open(file_name) as file:
        for line in file:
            if line.startswith('#'):
                continue
            else:
                Al.append(convert(line.split(sep=sep)[0]))
                Bl.append(convert(line.split(sep=sep)[1]))
    return Al, Bl


rA = write_list('V(H-H)_AK.txt')[0]#ovo je x
vK = write_list('V(H-H)_AK.txt')[1]#ovo je y

rA_dots = np.arange(rA[0], rA[-1], 0.1)

cs = CubicSpline(rA, vK)

cs(rA_dots)
