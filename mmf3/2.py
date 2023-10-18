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

filex1 = open("2x1.txt","a")
filex1.truncate(0)
cols=("h", "error")
filex1.write('\t'.join(cols)+'\n')
for l in range (1,7):
    h=10**(-l)
    cols=(str(h),str('{:e}'.format(error(d2fdx2_e_to_the_x(1,h),math.exp(1)))))
    filex1.write('\t'.join(cols) + '\n')

filex5 = open("2x5.txt","a")
filex5.truncate(0)
cols=("h", "error")
filex5.write('\t'.join(cols)+'\n')
for l in range (1,7):
    h=10**(-l)
    cols=(str(h),str('{:e}'.format(error(d2fdx2_e_to_the_x(5,h),math.exp(5)))))
    filex5.write('\t'.join(cols) + '\n')

filex10 = open("2x10.txt","a")
filex10.truncate(0)
cols=("h", "error")
filex10.write('\t'.join(cols)+'\n')
for l in range (1,7):
    h=10**(-l)
    cols=(str(h),str('{:e}'.format(error(d2fdx2_e_to_the_x(10,h),math.exp(10)))))
    filex10.write('\t'.join(cols) + '\n')
        