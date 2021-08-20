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
from datetime import datetime
import matplotlib.pyplot as plt
import sys
import csv

parser = argparse.ArgumentParser(description='A cool program.')
parser.add_argument('--dwave', dest='dwave', action='store_const', const=True, default=False, help='Flag to run on the dwave quantum annealer')

args = parser.parse_args()


def printLine():
    lineWidth = 60
    for i in range(lineWidth):
        print("-", end="")
    print("")
def feed(string, flush=True):
    lineWidth = 60
    remaining = lineWidth - len(string)
    pre = ""
    for i in range(math.floor(remaining/2)):
        pre += "-"
    post = ""
    for i in range(math.ceil(remaining/2)):
        post += "-"
    print(pre+string+post, end="\r")
    if(flush):
        sys.stdout.flush()
def createVariable(label, bits, scale):
    binary = sympy.symbols("%s:%s"%(label,bits))
    expr = "0"
    for i in range(bits-1):
        expr += "+ %s * 2**%s*%s"%(binary[i], i, scale)
    expr += "- %s * 2**%s*%s"%(binary[bits-1], bits-1, scale)
    expr = parse_expr(expr)
    return binary, expr
def get_variables(expr, constants):
    x = []
    for e in expr.free_symbols:
        if(e not in constants): #Don't count constants here
            x.append(e)
    return x
def reduceExpr(expr, constants): #Takes exponentiated binary terms and reduces them
    x = get_variables(expr, constants)
    toSubs = []
    degree = max(sympy.degree_list(expr))
    for var in x:
        for i in reversed(range(2,degree+1)):
            toSubs.append((var**i, var))
    return expr.subs(toSubs)
def getTime():
    #Purely for logging purposes
    return datetime.now().strftime("%H:%M")
def full_expand(expr, constants):
    expr = sympy.expand(expr)
    #Expression is the base expression to be expanded and quadritized. Example: expr=(N-pq)**2
    expr = reduceExpr(expr, constants)
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
        #expr = reduceExpr(expr, constants) #Shouldn't be necessary?
        expr_args = sympy.Add.make_args(expr)
        for arg in expr_args:
            x = get_variables(arg, constants)
            if(len(x) > 2): #Method from Pseudo-Boolean Optimization (Boros, Hammer)
                feed("%s: AncIndex: %s, x: %s"%(getTime(),ancIndex,x))
                anc = sympy.symbols(letter_index(ancIndex, "anc"))
                expr += M*x[0]*x[1]
                expr += -2*M*x[0]*anc
                expr += -2*M*x[1]*anc
                expr += 3*M*anc
                coeffs = expr.coeff(x[0]).coeff(x[1]) #Not efficient, but coeff(x[0]*x[1]) doesn't work
                toSubs = []
                for coeff in coeffs.args:
                    if(len(get_variables(coeff, constants)) > 0):
                        toSubs.append((coeff*x[0]*x[1],coeff*anc))
                expr = expr.subs(toSubs)
                ancIndex += 1
                quadritized = False
                if(len(x) > 3): #Means it needs multiple passes
                    expr = sympy.expand(expr)
                break
    feed("AncIndex: %s"%(ancIndex), flush=False)
    return expr
def storeForwards(bits, scale, label="basic"):
    f_cur, F_cur = createVariable("f_cur", bits, scale)
    x_cur, X_cur = createVariable("x_cur", bits, scale)
    x_prev, X_prev = createVariable("x_prev", bits, scale)
    f_other_prev, F_other_prev = createVariable("f_other_prev", bits, scale)
    const = {
        "gamma": sympy.symbols("gamma"),
        "Lambda": sympy.symbols("Lambda"),
        "beta": sympy.symbols("beta"),
        "C": sympy.symbols("C"),
    }
    var = {
        "f_cur": F_cur,
        "x_cur": X_cur,
        "x_prev": X_prev,
        "f_other_prev": F_other_prev,
    }
    binary = {
        "f_cur": f_cur,
        "x_cur": x_cur,
        "x_prev": x_prev,
        "f_other_prev": f_other_prev
    }
    #Introducing the sigmoid part (order doesn't matter)
    expr = (var["x_cur"]-(const["C"]+const["Lambda"]*var["x_prev"]-const["beta"]*var["f_other_prev"]))**2
    expr += ((1+(const["gamma"]*var["x_cur"])**2)*(var["f_cur"]-0.5)-const["gamma"]*var["x_cur"])**2
    constants = [const[key] for key in const.keys()]
    expr = full_expand(expr, constants)
    with open("pickles/%s[%s][%s][f].pkl"%(label,bits,scale), "wb") as f:
        pickle.dump([bits,scale,expr,binary,const,"f"],f)

def loadPolynomial(params, fileName):
    #Going to need to substitute in each of the timesteps properly. Currently labeled with cur for t, and prev for t-1
    #Also need to fill in the problem type. Currently is unlabeled for the default, and other for the other action
    with open(fileName, 'rb') as f:
        bits,scale,base_expr,binary,const,direction = pickle.load(f)
    #Creating a random set of points according to y=mx+b for random, integer m and b
    I = {'col': [], 'dir': []}
    C = {'col': [], 'dir': []} #Technically can figure one from the other, but I already created two columns in the csv
    actions = ['col', 'dir']
    if(direction == "f"):
        with open('in.csv') as csv_in:
            reader = csv.reader(csv_in)
            for row in reader:
                I['col'].append(row[0])
                I['dir'].append(row[1])
                C['col'].append(row[2])
                C['dir'].append(row[3])
    f = {'col': [], 'dir': []}
    x = {'col': [], 'dir': []}
    pts_used = 5 #Should be len(C['col']) for the full dataset
    expr = 0
    for i in range(pts_used):
        for action in actions:
            feed("Substituting '%s' at %s"%(action,i))
            f[action].append(sympy.symbols("%s:%s"%(letter_index(i,"f_%s"%(action)), bits)))
            x[action].append(sympy.symbols("%s:%s"%(letter_index(i,"x_%s"%(action)), bits)))
            toSubs = []
            for key in const.keys(): #Each key is the name of a constant
                if(key in params):
                    toSubs.append((const[key], params[key]))
                if(key == "C"):
                    toSubs.append((const[key], C[action][i]))
            for b in range(bits):
                toSubs.append((binary['f_cur'][b], f[action][i][b]))
                toSubs.append((binary['x_cur'][b], x[action][i][b]))
                toSubs.append((binary['x_prev'][b], 0 if i==0 else x[action][i-1][b]))
                toSubs.append((binary['f_other_prev'][b], 0 if i==0 else f['dir' if action=='col' else 'col'][i-1][b]))
            expr += base_expr.subs(toSubs)
    feed("Substitution complete!", flush=False)

    bqm = dimod.BinaryQuadraticModel({}, {}, 0, dimod.BINARY)
    for var in expr.free_symbols: #Adds all the variables
        bqm.add_variable(str(var))
    for arg in expr.args:
        free = [e for e in arg.free_symbols]
        if(len(free) == 2):
            bqm.set_quadratic(str(free[0]),str(free[1]), float(arg.coeff(free[0]*free[1])))
        if(len(free) == 1):
            bqm.set_linear(str(free[0]), float(arg.coeff(free[0])))
        if(arg.is_constant()): #Handling the one constant shift
            bqm.add_offset(int(arg))
    print("Num. Variables: %s"%(len(bqm.variables)))

    returnDict = anneal(bqm, dwave=args.dwave, num_reads=1000, num_sweeps=1000000, cores=7, best=True)
    sampleset = list(returnDict.values())[0]
    sample = sampleset.first.sample
    print("Energy: %s"%(sampleset.first.energy))
    def convertBinary(sample, labels):
        total = 0
        num_bits = len(labels)
        for i in range(num_bits-1):
            total += sample[str(labels[i])]*2**(i)*scale
        total -= sample[str(labels[num_bits-1])]*2**(num_bits-1)*scale
        digits_after_decimal = abs(round(math.log(scale,10)))
        total = round(total, digits_after_decimal)
        return total
    f_result = {'col': [], 'dir': []}
    for i in range(pts_used):
        for action in actions:
            f_result[action].append(convertBinary(sample, f[action][i]))
    x_result = {'col': [], 'dir': []}
    for i in range(pts_used):
        for action in actions:
            x_result[action].append(convertBinary(sample, x[action][i]))
    #print(sample)
    print("f: %s"%(f_result))
    print("x: %s"%(x_result))
def main():
    startTime = time.time()
    bits = 6
    scale = 0.1
    #storeForwards(bits, scale)
    params = {'omega': 0.001, 'g': 1, 'lambda': 2, 'beta': 1, 'w': 0}
    params['gamma'] = params['g']/1
    params['Lambda'] = 1-params['lambda']
    loadPolynomial(params, "pickles/%s[%s][%s][%s].pkl"%("basic",bits,scale,"f"))
    printLine()
    print("Runtime (s): %s"%(time.time()-startTime))
if(__name__ == "__main__"):
    main()
