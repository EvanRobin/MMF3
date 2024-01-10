import math
import numpy as np

yo1 = 5
yo2 = 0.325
A = 1
B = 3
C = 2
D = 0.5

def y1(t):
    return yo1 + A*math.cos(B*t)

def y2(t):
    return yo2 + C*math.exp(D*t)
def y3(t):
    return math.cos(t)
#establish an interval a and b
#the function will take in the interval and the tolerance
def y3bisection(a,b,tolerance):
    ta = a
    tb = b
    if y1(ta)*y1(tb) < 0:
        while (np.abs(ta-tb)>=tolerance):
            tc = (ta+tb)/2 #c is the interval which will always be right between a and b
            #to check if the interval even changes signs 
            if y1(ta)*y1(tc) > 0: 
                #this means that a to c doesnt cointain the zero of the function therefore we should look from c to b
                ta = tc
            elif y1(tb)*y1(tc) > 0:
                #in this case the interval should be a to c
                tb = tc
            else:
                print("ta: ", ta, "tb: ", tb, "tc ", tc)
                ta = 0
                tb = 0
    else:
        ta = 0
        tb = 0
        print("Pick another interval")
        tc=0
    #if were out of the while loop then either something has gone wrong or it worked
    return tc

def y3bisection(a,b,tolerance):
    ta = a
    tb = b
    if y3(ta)*y3(tb) < 0:
        while (np.abs(ta-tb)>=tolerance):
            tc = (ta+tb)/2 
            if y3(ta)*y3(tc) > 0: 
                ta = tc
            elif y3(tb)*y3(tc) > 0:
                tb = tc
            else:
                print("ta: ", ta, "tb: ", tb, "tc ", tc)
                ta = 0
                tb = 0
    else:
        ta = 0
        tb = 0
        print("Pick another interval")
        tc=0
    return tc

#for the newton rhapsody method the derivitive canot be 0

def newtrhapsmethod(func, funcderivetive, x, n):

    def f(x):
        f = eval(func)#eval is a built-in- function used in python, eval function parses the expression argument and evaluates it as a python expression. In simple words, the eval function evaluates the “String” like a python expression and returns the result as an integer
        return f
    def df(x):
        df = eval(funcderivetive)#OVO JE PRE DOBRO
        return df
    
    for i in range(1,n):#how many times we search for a closer aproxiamtion
        xo = x - (f(x)/df(x))
        x=xo #we are updating x for the next aproximation
    return x
'''print(newtrhapsmethod("x**2 - 2","x*2", 10, 100))'''

def generalbisection(func,a,b,tolerance):
    def f(x):
        f = eval(func)
        return f
    ta = a
    tb = b
    if f(ta)*f(tb) < 0:
        while (np.abs(ta-tb)>=tolerance):
            tc = (ta+tb)/2 
            if f(ta)*f(tc) > 0: 
                ta = tc
            elif f(tb)*f(tc) > 0:
                tb = tc
            else:
                print("ta: ", ta, "tb: ", tb, "tc ", tc)
                ta = 0
                tb = 0
    else:
        ta = 0
        tb = 0
        print("Pick another interval")
        tc=0
    return tc

print(generalbisection("math.cos(x)",1,2,0.004))
