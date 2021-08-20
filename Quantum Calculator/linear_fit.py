import argparse
from bqm_calculator import BQMCalculator
from multicore_anneal import anneal
from twos_complement import toDecimal, toBinary
import matplotlib.pyplot as plt
from random import randint
import math

bits = 6

parser = argparse.ArgumentParser(description='A cool program.')
parser.add_argument('--dwave', dest='dwave', action='store_const', const=True, default=False, help='Flag to run on the dwave quantum annealer')
args = parser.parse_args()

#Creating a random set of points according to y=mx+b for random, integer m and b
max_num = 2**(bits-1)-1
m = 1
b = 3
x = []
max_x = math.floor((max_num-b)/m)
for i in range(2):
    r = randint(-max_x, max_x)
    if(r not in x):
        x.append(r)
y = [m*X+b for X in x]
plt.plot(x, y, 'bo')
plt.xlabel("X")
plt.ylabel("Y")
plt.xlim(-max_num, max_num)
plt.ylim(-max_num, max_num)
plt.savefig("out.png")

calculator = BQMCalculator(bits)

calculator.addNum("m", m)
calculator.addNum("b", b)
for i in range(len(x)):
    calculator.addNum("x{%s}"%(i), x[i])
    calculator.addNum("y{%s}"%(i), y[i])
    calculator.multiply("x{%s}"%(i), "m", "mx{%s}"%(i))
    calculator.add("mx{%s}"%(i), "b", "y{%s}"%(i))
bqm = calculator.getBQM()

def printSample(sampleset, label):
    lineWidth = 40
    for i in range(lineWidth):
        print("-", end="")
    print("-")
    for i in range(int((lineWidth-len(label))/2)):
        print(" ", end="")
    print("[%s]"%(label))
    #print(sample)
    print("Energy: %s"%(sampleset.first.energy))
    calculator.report(sampleset.first.sample)
returnDict = anneal(bqm, dwave=args.dwave, num_reads=1000, num_sweeps=1000, best=True, cores=8)
for label in returnDict.keys():
    printSample(returnDict[label], label)
sampleset = list(returnDict.values())[0]
sample = sampleset.first.sample

m_predict = toDecimal(sample, "m", bits)
b_predict = toDecimal(sample, "b", bits)
#m_predict = m
#b_predict = b
plt.plot([-max_num, max_num], [-m_predict*max_num+b_predict, m_predict*max_num+b_predict])
plt.savefig("out.png")

print("Energy: %s"%(sampleset.first.energy))
print("Predicted: y={}x+{}".format(m_predict, b_predict))
print("Actual: y={}x+{}".format(m, b))
