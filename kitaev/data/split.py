from numpy import loadtxt as ltxt
import sys

try:
	str_input = sys.argv[1]
except OSError:
	print('Error: input file')
	raise
#op = open(file = str_input, mode='r')

x = ltxt(str_input, unpack = 'True')

str_out1 = '%s1' %(str_input)
str_out2 = '%s2' %(str_input)

try:
	f = open(file=str_out1, mode='w')
	g = open(file=str_out2, mode='w')
except:
	f = open(file=str_out1, mode='x')
	g = open(file=str_out2, mode='x')

for j in range(0, len(x[:,0])):
	f.write('%.16e	' %(x[j,0]) )
f.write('\n')

for i in range(0, len(x[0,:])-1):
	if(x[0, i+1] < 0.):
		for j in range(0, len(x[:, i+1])):
			g.write('%.16e	' %(x[j, i+1]) )
		g.write('\n')
	else:
		for j in range(0, len(x[:, i+1])):	
			f.write('%.16e	' %(x[j, i+1]) )
		f.write('\n')

f.close()
g.close()
'''
from numpy import loadtxt as ltxt
import sys

try:
	str_input = sys.argv[1]
except OSError:
	print('Error: input file')
	raise
#op = open(file = str_input, mode='r')

try:
	x, y, w, z = ltxt(str_input, unpack = 'True')
except A:
	x, y, w = ltxt(str_input, unpack = 'True')
except B:
	x, y, w, s, z = ltxt(str_input, unpack = 'True')

str_out1 = '%s1' %(str_input)
str_out2 = '%s2' %(str_input)

try:
	f = open(file=str_out1, mode='w')
	g = open(file=str_out2, mode='w')
except:
	f = open(file=str_out1, mode='x')
	g = open(file=str_out2, mode='x')
try:
	f.write('%.16e	%.16e	%.16e	%.16e\n' %(x[0], y[0], w[0], z[0]) )
except A:
		f.write('%.16e	%.16e	%.16e\n' %(x[0], y[0], w[0]) )
except B:
	f.write('%.16e	%.16e	%.16e	%.16e	%.16e\n' %(x[0], y[0], w[0], s[0], z[0]) )


for i in range(0, len(x)-1):
	if(x[i+1] < x[i]):
		try:
			g.write('%.16e	%.16e	%.16e	%.16e\n' %(x[i+1], y[i+1], w[i+1], z[i+1]) )
		except A:
				g.write('%.16e	%.16e	%.16e\n' %(x[i+1], y[i+1], w[i+1]) )
		except B:
			g.write('%.16e	%.16e	%.16e	%.16e	%.16e\n' %(x[i+1], y[i+1], w[i+1], s[i+1], z[i+1]) )
	else:
		try:
			f.write('%.16e	%.16e	%.16e	%.16e\n' %(x[i+1], y[i+1], w[i+1], z[i+1]) )
		except A:
				f.write('%.16e	%.16e	%.16e\n' %(x[i+1], y[i+1], w[i+1]) )
		except B:
			f.write('%.16e	%.16e	%.16e	%.16e	%.16e\n' %(x[i+1], y[i+1], w[i+1], s[i+1], z[i+1]) )

f.close()
g.close()
'''
