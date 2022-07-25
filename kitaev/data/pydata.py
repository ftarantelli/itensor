from numpy import *


x,y,z,w = loadtxt('wireABCu10k1l400.dat', unpack = 'True')

diff = abs( w[0] - w[1])/2.

for i in range(0,len(x)):
	if abs(y[i]-w[i]) < diff :
		print("Omega = ", x[i], y[i], w[i], "A = ", z[i] )
        
for i in range(0,len(x)):
	if abs(y[i]-w[i]) < diff :
		print('set arrow from %.10f, graph 0 to %.10f,graph 1 nohead lt 8 dt 2 lw 3\n' %(x[i]))
