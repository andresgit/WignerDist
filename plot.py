import matplotlib.pyplot as plt
import math
import numpy

file = open("test1.txt")
interval = float(file.readline())
data = [list(map(float,x.strip("()").split(","))) for x in file.read().split()]
file.close()
##print(data)

print(interval)

N = int(math.sqrt(len(data)))

##mag = [x[0]*x[0]+x[1]*x[1] for x in data]
mag = [x[0] for x in data]
##print(mag)
Z = numpy.array(mag).reshape((N,N))
print(Z)

x = numpy.linspace(-interval/2, interval/2,N)
p = numpy.linspace(-math.pi*N/(2*interval),math.pi*N/(2*interval),N)
X,P = numpy.meshgrid(x,p,indexing='ij')

##print(X)
##print(P)

plt.figure()
plt.contourf(X, P, Z)
plt.ylim(-interval/2,interval/2)
plt.colorbar()
plt.savefig("test1.png")
plt.show()

