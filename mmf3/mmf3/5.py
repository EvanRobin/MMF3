import matplotlib.pyplot as plt
import numpy as np
import math

x = [1, 2, 3, 4, 5]
y = [-5, -12, -15, -8, 15]

def l(x, y, t):
    p = 0
    for i in range(0,len(x)):
        l=1
        for j in range(0,len(x)):
            if j != i:
                l *= (t-x[j])/(x[i]-x[j])
        p += l*y[i]
    return p
'''
t = np.arange(1, 5, 0.001)
p = []
for i in range (0,len(t)):
    p.append(l(x, y, t[i]))
'''
'''  
plt.plot(t, p)  
plt.xlabel('x - axis')  
plt.ylabel('y - axis')  
plt.title('My first graph!')  
plt.show()  
''' 

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


rA = write_list('V(H-H)_AK.txt')[0]
vK = write_list('V(H-H)_AK.txt')[1]

rA_dots = np.arange(rA[0], rA[-1], 0.1)
pvK = []
for i in range (0,len(rA_dots)):
    pvK.append(l(rA, vK, rA_dots[i]))



import sys
# polint radi polinomnu interpolaciju na temelju Nevilleva algoritma
# ulaz: NP parova podataka (xi,yi) i argument x
# izlaz: vrijednost polinoma yN u tocki x i greska procjene dy
def polint(xi, yi, NP, x):
	ns=1; C=[0]; D=[0]; xa=[0]; ya=[0];
	for i in range(NP):
		xa.append(xi[i])
		ya.append(yi[i])
	# trazimo najblizeg susjeda (ns) od x
	mdx=abs(x-xa[1]) # udaljenost tocke x od 1. cvora
	for i in range(1,NP+1):
		dx=abs(x-xa[i])
		if (dx < mdx): # udaljenost od ostalih cvorova
			ns=i
			mdx=dx # minimala udaljenost
		# pocetne vrijednosti (nulti stupac)
		C.append(ya[i])
		D.append(ya[i])
	yN=ya[ns] # prva aproksimacija (prvi stupac)
	ns=ns-1
	for m in range(1,NP):         # za svaki stupac
		for i in range(1,NP-m+1): # za svaki redak
			bCx=xa[i]-x         # brojnik u C: razlika x-eva
			bDx=xa[i+m]-x       # brojnik u D: razlika x-eva
			CD=C[i+1]-D[i]       # razlika iz prethodnog stupca
			# stop ako postoje zaokruzeno-isti xi
			odnos=bCx-bDx
			if (odnos == 0.0): sys.exit("STOP: dijeljenje s 0!")
			odnos=CD/odnos
			D[i]=bDx*odnos  # koeficijenti C i D
			C[i]=bCx*odnos
		if (2*ns < (NP-m)):
			dy=C[ns+1]
		else:
			dy = D[ns]
			ns = ns - 1
		yN += dy
	return yN, dy
#yN je vrijednost interpol polinoma a dy je greska

pvKP = []
error = []
for i in range (0,len(rA_dots)):
    pvKP.append(polint(rA, vK, len(rA) ,rA_dots[i])[0])
    error.append(abs(polint(rA, vK, len(rA) ,rA_dots[i])[1]))
print(error)
plt.scatter(rA_dots, pvK)
plt.scatter(rA_dots, pvKP) 
plt.errorbar(rA_dots, pvKP, error, fmt = 'o')
plt.ylim(-10, 10) 
plt.xlim(1, 10)
plt.xlabel('r/A')  
plt.ylabel('V/K')  
plt.title('<(")')  
plt.show() 