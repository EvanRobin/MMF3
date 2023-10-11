import math
import numpy as np

def e_to_the_x(x):
    return math.exp(x)

def d2fdx2_e_to_the_x(x,h):\
    return (e_to_the_x(x+h)+e_to_the_x(x-h)-2*e_to_the_x(x))/h**2
