import matplotlib.pyplot as plt
import math
import numpy

file = open("test1.txt")
interval = float(file.readline())
data = [float(x) for x in file.read().split()]
file.close()
##print(data)

N = int(math.sqrt(len(data)))

print("On graph, X points: {}, P points: {:.0f}".format(N, interval*interval/math.pi))

Z = numpy.array(data).reshape((N,N))
##print(Z)

x = numpy.linspace(-interval/2, interval/2,N)
p = numpy.linspace(-math.pi*N/(2*interval),math.pi*N/(2*interval),N)
X,P = numpy.meshgrid(x,p,indexing='ij')

##print(X)
##print(P)

plt.figure()
plt.contourf(X, P, Z)
plt.xlim(-interval/2,interval/2)
plt.ylim(-interval/2,interval/2)
plt.colorbar()
plt.savefig("test1.png")
plt.show()

