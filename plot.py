import matplotlib.pyplot as plt
import math
import numpy

def plot1D(filename):
    file = open("build/pongrid.txt")
    interval = float(file.readline())
    data = [complex(x.strip("()").replace(",","+").replace("+-","-")+"j") for x in file.read().split()]
    N = len(data)
    p1 = numpy.linspace(-math.pi*N/(interval),math.pi*N/(interval),N)

##    print(p1)

    plt.xlim(-3,3)
    datap = list(map(lambda x: abs(x)*abs(x),data))
    plt.plot(p1, datap, label="DFT")

    file = open("build/checkProb.txt")
    interval = float(file.readline())
    data = [float(x) for x in file.read().split()]
    N = int(math.sqrt(len(data)))
    x = numpy.linspace(-interval/2, interval/2,N)
    p = numpy.linspace(-math.pi*N/(2*interval),math.pi*N/(2*interval),N)
    Z = numpy.array(data).reshape((N,N))
    X,P = numpy.meshgrid(x,p,indexing='ij')

    xval = math.pi/interval*Z.sum(1)
    pval = interval/N*Z.sum(0)

##    print(x)

##    for i in range(N):
##        print(i,xval[i], 1/math.sqrt(math.pi)*numpy.exp(-x[i]*x[i]))
##    print()
    for i in range(N):
        print("{} {:.7f}  {:.7f}  {:.7f}  {:.7f}".format(i, p1[i], pval[i], datap[i], 1/math.sqrt(math.pi)*numpy.exp(-p1[i]*p1[i])))

    shift = 0.0
    plt.plot(x, xval, label="Wig x")
    plt.plot(p,pval+shift, label="Wig p")
    plt.plot(x, 1/math.sqrt(math.pi)*numpy.exp(-x*x)+2*shift, label="From wavefunc")
    plt.legend()
    plt.savefig(filename.replace(".txt",".png"))
    plt.show()
    

def plot(filename, show=True):
    file = open("build/{}".format(filename))
    interval = float(file.readline())
    data = [float(x) for x in file.read().split()]
    file.close()
    ##print(data)

    N = int(math.sqrt(len(data)))

    print("On graph, X points: {}, P points: {:.0f}".format(N, interval*interval/math.pi))

    Z = numpy.array(data).reshape((N,N))
##    print(Z)

    

    x = numpy.linspace(-interval/2, interval/2,N)
    p = numpy.linspace(-math.pi*N/(2*interval),math.pi*N/(2*interval),N)
    X,P = numpy.meshgrid(x,p,indexing='ij')

##    print(X)
##    print(P)

    xval = math.pi/interval*Z.sum(1)
    pval = interval/N*Z.sum(0)

    print(x)

    for i in range(N):
        print(i,xval[i], 1/math.sqrt(math.pi)*numpy.exp(-x[i]*x[i]))
    print()
    for i in range(N):
        print(i,pval[i], 1/math.sqrt(math.pi)*numpy.exp(-p[i]*p[i]))

    shift = 0.0
    plt.plot(x, xval)
    plt.plot(p,pval+shift)
    plt.plot(x, 1/math.sqrt(math.pi)*numpy.exp(-x*x)+2*shift)
    plt.show()
    
##    plt.figure()
##    plt.contourf(X, P, Z)
##    plt.xlim(-interval/2,interval/2)
##    plt.ylim(-interval/2,interval/2)
##    plt.colorbar()
##    plt.savefig(filename.replace(".txt",".png"))
##    if show:
##        plt.show()

filenames = [
##    "test1.txt",
##    "oscEigen0.txt",
##    "oscEigen1.txt",
##    "oscEigen2.txt",
##    "oscCoh00.txt",
##    "oscCoh10.txt",
##    "oscCoh01.txt",
##    "oscCoh11.txt",
##    "oscPure.txt",
##    "oscMixed.txt",
    "checkProb.txt",
    ]

plot1D("pongrid.txt")
##for filename in filenames:
##    plot(filename, len(filenames)==1)
