import matplotlib.pyplot as plt
import math
import numpy

file = open("test1.txt")
data = [list(map(float,x.strip("()").split(","))) for x in file.read().split()]
file.close()
print(data)


N = int(math.sqrt(len(data)))

mag = [x[0]*x[0]+x[1]*x[1] for x in data]
print(mag)
Z = numpy.array(mag).reshape((N,N))

x = numpy.linspace(0,1,N)
X,Y = numpy.meshgrid(x,x)
plt.figure()
plt.contour(Y, X, Z)
plt.show()

