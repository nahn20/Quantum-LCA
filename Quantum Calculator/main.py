import argparse
from bqm_calculator import BQMCalculator
from multicore_anneal import anneal

parser = argparse.ArgumentParser(description='A quantum binary calculator using Boolean logic')
parser.add_argument('--dwave', dest='dwave', action='store_const', const=True, default=False, help='Flag to run on the dwave quantum annealer')
NUM_BITS = 5 #In two's complement, so max would be 2**(NUM_BITS-1)

args = parser.parse_args()

calculator = BQMCalculator(NUM_BITS)
calculator.addNum("m", -4)
calculator.multiply("P", "n", "q", scale=1)

bqm = calculator.getBQM()
print("Number of Variables: %s"%(len(bqm.variables)))

def printBits(sample, var):
    i = 0
    while("%s[%s]"%(var,i) in sample.keys()):
        i += 1
    i -= 1
    while("%s[%s]"%(var,i) in sample.keys()):
        #print("'%s[%s]': %s"%(var,i,sample["%s[%s]"%(var,i)]), end = ", ")
        print("%s"%(sample["%s[%s]"%(var,i)]), end = "")
        i -= 1
    print(" | '%s'"%(var))
def printSample(sampleset, label): #Prints a single sample given from multicore_anneal
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
returnDict = anneal(bqm, dwave=args.dwave, num_reads=100, num_sweeps=1000)
for label in returnDict.keys():
    printSample(returnDict[label], label)
