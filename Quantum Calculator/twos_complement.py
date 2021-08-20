import math
#Binary example:
def sign(num):
    if(num >= 0):
        return 1
    else:
        return -1
#Returns as a dictionary with same formatting as below. 
def toDecimal(dictionary, var, num_bits):
    total = 0
    total -= dictionary["%s[%s]"%(var,num_bits-1)]*2**(num_bits-1)
    for i in range(num_bits-1):
        total += dictionary["%s[%s]"%(var,i)]*2**i
    return total

#Takes inputs as a dictionary of var0-varn, and var.0-var.n. Periods denote a split to decimal.
def toBinary(num, var, num_bits):
    binary = {}
    if(num < 0):
        num += 2**(num_bits-1)
        binary["%s[%s]"%(var,num_bits-1)] = 1
    else:
        binary["%s[%s]"%(var,num_bits-1)] = 0
    def binarify(value, var, i):
        exp = 2**i
        bit = value // exp
        binary["%s[%s]"%(var,i)] = int(bit)
        if(i == 0):
            return
        binarify(value-bit*exp, var, i-1)
    binarify(num, var, num_bits-2)
    return binary
