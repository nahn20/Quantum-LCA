import dimod
import argparse
import sympy
from sympy.parsing.sympy_parser import parse_expr
from letter_index import letter_index
import time
import pickle
from multicore_anneal import anneal
from dwave.system import LeapHybridSampler
import math
from random import randint
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='y=mx+b linear fit problem using quantum annealing')
parser.add_argument('--dwave', dest='dwave', action='store_const', const=True, default=False, help='Flag to run on the dwave quantum annealer')

args = parser.parse_args()


def printLine(): #Output formatting
    lineWidth = 40
    for i in range(lineWidth):
        print("-", end="")
    print("")
#See Quantum LCA for a more optimized version of this function
def storePolynomial(fileName, bits, scale):
    #Expands the problem polynomial and stores it in the given fileName
    #Bits is the number of bits used for each variable. Includes negative representation, so each variable is -2^(bits-1) <= x <= 2^(bits-1)
    #Scale is the shift scale for the problem. For example, a scale of 0.1 would make a problem with 5 bits range from -1.6 <= x <= 1.5
    def createVariable(label):
        binary = sympy.symbols("%s:%s"%(label,bits))
        expr = "0"
        for i in range(bits-1):
            expr += "+ %s * 2**%s*%s"%(binary[i], i, scale)
        expr += "- %s * 2**%s*%s"%(binary[bits-1], bits-1, scale)
        expr = parse_expr(expr)
        return binary, expr
    #binary is the sympy collection of the binary variables. expr is the combined binary variables with proper scales so that they fill the range
    slope, Slope = createVariable("m")
    intercept, Intercept = createVariable("b")
    x_coord = sympy.symbols("x")
    y_coord = sympy.symbols("y")
    expr = sympy.expand((y_coord - (Slope*x_coord+Intercept))**2) #Our problem equation
    for var in expr.free_symbols: #Can safely do this because we're working in binary. 1^2 = 1 and 0^2 = 0
        if(var != x_coord and var != y_coord): #Don't count constants here
            expr = expr.subs(var**2, var)
    #SymPy variable names can only use letters, so we use this ancIndex combined with a conversion to alphabet index
    ancIndex = 0
    quadritized = False
    #Only need to calculate M once then it's used for the rest of the cycles
    M = 1 #See Pseudo-Boolean Optimization (Boros, Hammer)
    poly = sympy.Poly(expr)
    for coeff in poly.coeffs():
        M += abs(coeff)
    while(not quadritized):
        quadritized = True
        expr_args = sympy.Add.make_args(expr)
        for arg in expr_args:
            x = []
            for e in arg.free_symbols:
                if(e != x_coord and e != y_coord): #Don't count constants here
                    x.append(e)
            if(len(x) > 2): #Method from Pseudo-Boolean Optimization (Boros, Hammer)
                anc = sympy.symbols(letter_index(ancIndex, "anc"))
                expr += M*x[0]*x[1]
                expr += -2*M*x[0]*anc
                expr += -2*M*x[1]*anc
                expr += 3*M*anc
                coeffs = expr.coeff(x[0]*x[1])
                for coeff in coeffs.args:
                    if(len(coeff.free_symbols) > 0):
                        expr = expr.subs(coeff*x[0]*x[1],coeff*anc)
                ancIndex += 1
                quadritized = False
            if(len(x) > 3): #Means it needs multiple passes
                expr = sympy.expand(expr)
    with open(fileName, 'wb') as f:
        pickle.dump([bits,scale,expr,y_coord,x_coord,slope,intercept],f)

def loadPolynomial(fileName):
    with open(fileName, 'rb') as f:
        bits,scale,base_expr,y_coord,x_coord,slope,intercept = pickle.load(f)
    #Creating a random set of points according to y=mx+b for random, integer m and b
    max_num = math.floor((2**(bits-1)-1)*scale)
    m = -0.5
    b = 200
    x = []
    for i in range(8):
        r = randint(-max_num, max_num)
        if(r not in x):
            x.append(r)
    y = [m*X+b for X in x]
    plt.plot(x, y, 'bo')
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.xlim(-max_num, max_num)
    plt.ylim(-max_num, max_num)
    plt.savefig("out.png")

    expr = 0
    for i in range(len(y)):
        #Each group represents the storePolynomial pickle for one set of variables. We need on of those for each of the points
        group = base_expr.subs(y_coord, y[i])
        group = group.subs(x_coord, x[i])
        expr += group
    bqm = dimod.BinaryQuadraticModel({}, {}, 0, dimod.BINARY)
    for var in expr.free_symbols: #Adds all the variables
        bqm.add_variable(str(var))
    for arg in expr.args:
        x = [e for e in arg.free_symbols]
        if(len(x) == 2):
            bqm.set_quadratic(str(x[0]),str(x[1]), float(arg.coeff(x[0]*x[1])))
        if(len(x) == 1):
            bqm.set_linear(str(x[0]), float(arg.coeff(x[0])))
        if(arg.is_constant()): #Handling the one constant shift
            bqm.add_offset(int(arg))
    print("Num. Variables: %s"%(len(bqm.variables)))

    returnDict = anneal(bqm, dwave=args.dwave, num_reads=10, num_sweeps=100, cores=6, best=True)
    sampleset = list(returnDict.values())[0]
    sample = sampleset.first.sample
    print("Energy: %s"%(sampleset.first.energy))
    def convertBinary(sample, labels):
        total = 0
        num_bits = len(labels)
        for i in range(num_bits-1):
            total += sample[str(labels[i])]*2**(i)*scale
        total -= sample[str(labels[num_bits-1])]*2**(num_bits-1)*scale
        return total
    slope_predict = convertBinary(sample, slope)
    intercept_predict = convertBinary(sample, intercept)
    print("Predict: Y=%s*X+%s"%(slope_predict, intercept_predict))
    print("Actual: Y=%s*X+%s"%(m, b))
    plt.plot([-max_num, max_num], [-slope_predict*max_num+intercept_predict, slope_predict*max_num+intercept_predict])
    plt.savefig("out.png")

def main():
    startTime = time.time()
    bits = 16
    scale = 0.01
    #storePolynomial("pickles/linear_%s_%s.pkl"%(bits,scale), bits, scale)
    loadPolynomial("pickles/linear_%s_%s.pkl"%(bits,scale))
    printLine()
    print("Runtime (s): %s"%(time.time()-startTime))
if(__name__ == "__main__"):
    main()
