import math
import sys
import matplotlib.pyplot as plt

filename = sys.argv[1]
outFileName = sys.argv[2]


class MSMCresult:
    """
    Simple class to read a MSMC result file. Constructor takes filename of MSMC result file
    """
    def __init__(self, filename):
        f = open(filename, "r")
        self.times_left = []
        self.times_right = []
        self.lambdas = []
        f.next() # skip header
        for line in f:
            fields = line.strip().split()
            nr_lambdas = len(fields) - 3
            if len(self.lambdas) == 0:
                self.lambdas = [[] for i in range(nr_lambdas)]
            time_left = float(fields[1])
            time_right = float(fields[2])
            self.times_left.append(time_left)
            self.times_right.append(time_right)
            for i in range(nr_lambdas):
                l = float(fields[3 + i])
                self.lambdas[i].append(l)
        self.T = len(self.times_left)

def popSizeStepPlot(filename, mu=1.25e-8, gen=30):
    """
    to be used with a step-plot function, e.g. matplotlib.pyplot.steps.
    returns (x, y), where x contains the left point of each step-segment in years, and y contains the effective population size. Note that there are two ways to make a step-plot. You should make sure that your step-plot routine moves right and then up/down instead of the other way around.
    If plotted on a logarithmic x-axis, you should adjust x[0] = x[1] / 4.0, otherwise the leftmost segement will start at 0 and won't be plotted on a log-axis.
    
    Options:
        mu: Mutation rate per generation per basepair (default=1.25e-8)
        gen: generation time in years (default=30)
    """
    M = MSMCresult(filename)
    x = [t * gen / mu for t in M.times_left]
    y = [(1.0 / l) / (2.0 * mu) for l in M.lambdas[0]]
    return (x, y)   
 
coordinate = popSizeStepPlot(filename, mu=3e-8, gen=1)

outFile = open(outFileName, 'w')

print coordinate
for num in range(0, len(coordinate[0])):
    outFile.write(str(coordinate[0][num]))
    outFile.write("\t")
    outFile.write(str(coordinate[1][num])+"\n")
    
#plt.plot(coordinate[0], coordinate[1])
plt.plot(coordinate[0][6:], coordinate[1][6:])
#plt.savefig("Mex_Lowland.png", dpi="300")  
plt.show()
  
    
    

