import math

file1 = open("2.txt","a")
file1.truncate(0) #ovo brise
cols=("x", "10**-1", "error", "10**-2", "error", "10**-3", "error", "10**-4", "error", "10**-5", "error", "10**-6", "error")
file1.write('\t'.join(cols)+'\n')

def e_to_the_x(x):
    return math.exp(x)

def d2fdx2_e_to_the_x(x,h):\
    return (e_to_the_x(x+h)+e_to_the_x(x-h)-2*e_to_the_x(x))/h**2

def error(a,b):
    return abs(a-b)


for i in range (0,11):
    line=[i]
    for l in range (1,7):
        h=10**(-l)
        line.append('{:e}'.format(d2fdx2_e_to_the_x(i,h)))
        line.append('{:e}'.format(error(d2fdx2_e_to_the_x(i,h),math.exp(i))))
    file1.write('\t'.join(str(item) for item in line) + '\n')

e_na_1=[]
e_na_5=[]
e_na_10=[]
for l in range (1,7):
    h=10**(-l)
    e_na_1.append(error(d2fdx2_e_to_the_x(1,h),math.exp(i)))
    e_na_5.append(error(d2fdx2_e_to_the_x(5,h),math.exp(i)))
    e_na_10.append(error(d2fdx2_e_to_the_x(10,h),math.exp(i)))