import math

file1 = open("1.txt","a")
file1.truncate(0) #ovo brise
cols=("x", "+-red(A)", "rekurzija(B)", "1/+red(C)", "clanova")
file1.write('\t'.join(cols)+'\n')
def exp_series(x,epsilon):
    k=1
    result = 1.0
    member = 1.0
    while abs(member) > epsilon:
        member *= -x/k
        result += member
        k += 1
    return k, result

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

for i in range (0,110,10): 
    a,b=exp_series(i,eps)
    c,d=exp_recursion(i,eps)
    e,f=exp_to_the_power_of(i,eps)
    cols1=(str(i),str('{:e}'.format(b)),str('{:e}'.format(d)),str('{:e}'.format(f)),str(a))
    file1.write('\t'.join(cols1)+'\n')

