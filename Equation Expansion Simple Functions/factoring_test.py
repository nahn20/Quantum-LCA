import dimod
import argparse
import sympy
from sympy.parsing.sympy_parser import parse_expr
from letter_index import letter_index
import time
import pickle
from multicore_anneal import anneal
from dwave.system import LeapHybridSampler

parser = argparse.ArgumentParser(description='A cool program.')
parser.add_argument('--dwave', dest='dwave', action='store_const', const=True, default=False, help='Flag to run on the dwave quantum annealer')

args = parser.parse_args()

def printLine():
    lineWidth = 40
    for i in range(lineWidth):
        print("-", end="")
    print("")
def storePolynomial(fileName):
    #Variables are only used to set up the equation and print out solutions
    a = sympy.symbols('a:2') #Setting up appropriate number of bits
    b = sympy.symbols('b:2')
    c = sympy.symbols('c:2')
    #Setting up p and q to use the above defined binary
    p = "0"
    q = "0"
    g = "0"
    for i in range(len(a)):
        p += "+ %s * 2**%s"%(a[i], i)
    for i in range(len(b)):
        q += "+ %s * 2**%s"%(b[i], i)
    for i in range(len(c)):
        g += "+ %s * 2**%s"%(c[i], i)
    p = parse_expr(p)
    q = parse_expr(q)
    g = parse_expr(g)
    N = sympy.symbols("N")
    #N = 7
    expr = sympy.expand((N - (p*q+g))**2) #Our problem equation
    for var in expr.free_symbols: #Can safely do this because we're working in binary
        if(var != N): #Don't count constants here
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
                if(e != N): #Don't count constants here
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
        pickle.dump([expr,N,a,b,c],f)

def loadPolynomial(fileName):
    with open(fileName, 'rb') as f:
        expr,N,a,b,c = pickle.load(f)
    expr = expr.subs(N, 317)
    bqm = dimod.BinaryQuadraticModel({}, {}, 0, dimod.BINARY)
    for var in expr.free_symbols: #Adds all the variables
        bqm.add_variable(str(var))
    for arg in expr.args:
        x = [e for e in arg.free_symbols]
        if(len(x) == 2):
            bqm.set_quadratic(str(x[0]),str(x[1]), int(arg.coeff(x[0]*x[1])))
        if(len(x) == 1):
            bqm.set_linear(str(x[0]), int(arg.coeff(x[0])))
        if(arg.is_constant()): #Handling the one constant shift
            bqm.add_offset(int(arg))
    print("Num. Variables: %s"%(len(bqm.variables)))
    print(bqm)
    def printSample(sampleset, label):
        sample = sampleset.first.sample
        printLine()
        print("[%s] Energy: %s"%(label, sampleset.first.energy))
        pTotal = 0
        for i in range(len(a)):
            pTotal += sample[str(a[i])]*2**(i)
        qTotal = 0
        for i in range(len(b)):
            qTotal += sample[str(b[i])]*2**(i)
        gTotal = 0
        for i in range(len(c)):
            gTotal += sample[str(c[i])]*2**(i)
        print("[%s] (%s)(%s)+%s = %s"%(label, pTotal, qTotal, gTotal, 317))
    returnDict = anneal(bqm, dwave=args.dwave, num_reads=1000, num_sweeps=10000, cores=6)
    for label in returnDict.keys():
        printSample(returnDict[label], label)
def main():
    startTime = time.time()
    #storePolynomial("linear_2.pkl")
    loadPolynomial("linear_5.pkl")
    printLine()
    print("Runtime (s): %s"%(time.time()-startTime))
if(__name__ == "__main__"):
    main()
