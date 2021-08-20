import dimod
from twos_complement import toDecimal, toBinary
from termcolor import colored

class BQMCalculator:
    #Stores a BQM with all the couplings and biases created by various gates
    #Since they're couplings and biases, they don't have labels
    #Auxiliary variables follow "(var0)(var1)description"
    def __init__(self, num_bits):
        self.bqm = dimod.BinaryQuadraticModel({}, {}, 0, dimod.BINARY)
        self.num_bits = num_bits #Not necessarily a constant, changes for the multiplication bit doubling case
        self.NUM_BITS = num_bits #num_bits but doesn't change
        self.checks = [] #A list of equations to check later to ensure they were calculated properly
        self.bitReport = [] #A collection of bits to report the value of
        self.numReport = []
    def report(self, sample): #Reports on all the self.checks
        #Formatted ["+", "a", "b", "c"]
        #Formatted ["=", "a", 10]
        def tfFormat(boolean):
            if(boolean):
                return colored("True ", "green")
            return colored("False", "red")
        for check in self.checks:
            if(check[0] == "+"):
                values = [toDecimal(sample,check[1],self.NUM_BITS), toDecimal(sample,check[2],self.NUM_BITS), toDecimal(sample,check[3],self.NUM_BITS)]
                overflow = sample["(%s)(%s)add-overflow"%(check[1],check[2])]
                print("%s | %s + %s == %s"%(tfFormat(values[0]+values[1]==values[2]),check[1], check[2], check[3]))
                print("\t%s + %s = %s"%(values[0], values[1], values[2]))
                if(abs(values[0]+values[1]) >= 2**(self.NUM_BITS-1)):
                    print("\t%s + %s out of bit range (%s)"%(values[0],values[1],self.NUM_BITS-1))
            if(check[0] == "*"):
                values = [toDecimal(sample,check[1],self.NUM_BITS), toDecimal(sample,check[2],self.NUM_BITS), toDecimal(sample,check[3],self.NUM_BITS)]
                print("%s | %s * %s == %s"%(tfFormat(values[0]*values[1]==values[2]),check[1], check[2], check[3]))
                print("\t%s * %s = %s"%(values[0], values[1], values[2]))
                if(abs(values[0]*values[1]) >= 2**(self.NUM_BITS-1)):
                    print("\t%s * %s out of bit range (%s)"%(values[0],values[1],self.NUM_BITS-1))
            if(check[0] == "="):
                value = toDecimal(sample,check[1],self.NUM_BITS)
                print("%s | %s == %s"%(tfFormat(value==check[2]),check[1],check[2]))
            if(check[0] == "!="):
                value = toDecimal(sample,check[1],self.NUM_BITS)
                print("%s | %s != %s"%(tfFormat(value!=check[2]),check[1],check[2]))
        for bit in self.bitReport:
            print("%s [%s]"%(sample[bit], bit))
        for var in self.numReport:
            for i in reversed(range(self.NUM_BITS)):
                print(sample["%s[%s]"%(var, i)], end="")
            print(" [%s]"%(var))
    def reportBit(self, var):
        self.bitReport.append(var)
    def reportNum(self, var):
        self.numReport.append(var)
    def add(self, in0, in1, out, overflow=False, negativeFlipAround=True, scale=1):
        #if(not self.bqm.has_variable("%s[0]"%(in0)) or not self.bqm.has_variable("%s[0]"%(in1))):
        #    raise ValueError("Binary variables %s and %s not found. Use calculator.addNum(var, num)"%(in0,in1))
        #Ignoring Two's Complement overflow for now. Perhaps eventually fix
        if(negativeFlipAround): #This is the flag for being used in multiplication, making sure it doesn't get too spammy
            self.checks.append(["+", in0, in1, out])
        carry = "(%s)(%s)add-carry"%(in0,in1)
        signCombine = "(%s)(%s)add-signCombine"%(in0,in1) #Controls the carry into the sign bit
        signCarry = "(%s)(%s)add-signCarry"%(in0,in1) #Works with signCombine
        auxSum = "(%s)(%s)add-auxSum"%(in0,in1) #Original sum without sign correction
        loopAround = "(%s)(%s)add-loopAround"%(in0,in1) #Corrects the sign
        bigger = in0
        if(not self.bqm.has_variable("%s[0]"%(in0))):
            self._addBinary(in0, self.num_bits)
        if(not self.bqm.has_variable("%s[0]"%(in1))):
            self._addBinary(in1, self.num_bits)
        if(not self.bqm.has_variable("%s[0]"%(out))):
            self._addBinary(out, self.num_bits)
        self._addBinary(carry, self.num_bits)
        self._halfAdder("%s[0]"%(in0), "%s[0]"%(in1), "%s[0]"%(out), "%s[0]"%(carry), scale=1*scale)
        for i in range(1, self.num_bits-1):
            self._fullAdder("%s[%s]"%(in0,i), "%s[%s]"%(in1,i), "%s[%s]"%(carry,i-1), "%s[%s]"%(out,i), "%s[%s]"%(carry,i), scale=1*scale)
        self._or("%s[%s]"%(in0,self.num_bits-1), "%s[%s]"%(in1,self.num_bits-1), signCombine, scale=1*scale)
        self._and("%s[%s]"%(carry,self.num_bits-2),signCombine,signCarry,scale=1*scale)
        #self._fullAdder("%s[%s]"%(in0,self.num_bits-1), "%s[%s]"%(in1,self.num_bits-1), signCarry, auxSum, loopAround, scale=1*scale)
        if(negativeFlipAround):
            self._fullAdderNoCarryOut("%s[%s]"%(in0,self.num_bits-1), "%s[%s]"%(in1,self.num_bits-1), signCarry, auxSum, scale=1*scale)
            self._and("%s[%s]"%(in0,self.num_bits-1), "%s[%s]"%(in1,self.num_bits-1), loopAround, scale=1*scale)
            self._or(auxSum, loopAround, "%s[%s]"%(out,self.num_bits-1), scale=1*scale)
        else:
            self._fullAdderNoCarryOut("%s[%s]"%(in0,self.num_bits-1), "%s[%s]"%(in1,self.num_bits-1), signCarry, "%s[%s]"%(out,self.num_bits-1), scale=1*scale)
        #Handling overflow
        overflow = "(%s)(%s)add-overflow"%(in0,in1)
        maxBitOff = "(%s)(%s)add-maxBitOff"%(in0,in1)
        maxBitOn = "(%s)(%s)add-maxBitOn"%(in0,in1)
        signBitOff = "(%s)(%s)add-signBitOff"%(in0,in1)
        signBitOn = "(%s)(%s)add-signBitOn"%(in0,in1)
        overflowPos = "(%s)(%s)add-overflowPos"%(in0,in1)
        overflowNeg = "(%s)(%s)add-overflowNeg"%(in0,in1)
        self._bothOff("%s[%s]"%(in0,self.num_bits-2),"%s[%s]"%(in1,self.num_bits-2),maxBitOff,scale=1*scale)
        self._and("%s[%s]"%(in0,self.num_bits-1),"%s[%s]"%(in1,self.num_bits-1),signBitOn,scale=1*scale)
        self._and(maxBitOff,signBitOn,overflowNeg,scale=1*scale)
        self._bothOff("%s[%s]"%(in0,self.num_bits-1),"%s[%s]"%(in1,self.num_bits-1),signBitOff,scale=1*scale)
        self._and("%s[%s]"%(in0,self.num_bits-2),"%s[%s]"%(in1,self.num_bits-2),maxBitOn,scale=1*scale)
        self._and(maxBitOn,signBitOff,overflowPos,scale=1*scale)
        self._or(overflowNeg,overflowPos,overflow,scale=1*scale)
        #Not currently punishing overflow here
    def multiply(self, in0, in1, out, overflow=False, scale=1):
        #See the fixing overflow pdf for examples

        #Rows are read start at 0 on the right (same as qubits)
        #Rows are summed together using regular addition
        #Rows are labeled (in0)(in1)mult-row[n][qubit]
        #Intermediary outcomes are labeled (in0)(in1)mult-combined[+1][qubit} for combined rows up to 1 (so rows 0 and 1)
        self.checks.append(["*", in0, in1, out])
        #Checks to see if the variable exists are done within _extend
        if(not self.bqm.has_variable("%s[%s]"%(in0, self.num_bits*2-1))): #Checks so it's not double extended
            self._extend(in0, scale=10*scale) #Locks the first num_bits+1 digits to be the same, effectively extending it
        if(not self.bqm.has_variable("%s[%s]"%(in1, self.num_bits*2-1))): #Checks so it's not double extended
            self._extend(in1, scale=10*scale)
        self.reportNum(in1)
        self.reportNum(in0)
        self.num_bits *= 2 #To detect overflow correctly, we need to extend our bits to twice as many
        for i in range(self.num_bits): #Going to do a more thorough check than usual. Same as in _extend, but we don't want to lock the output's qubits (since it could have overflow)
            if(not self.bqm.has_variable("%s[%s]"%(out, i))):
                self.addVar("%s[%s]"%(out,i))
        row = ["(%s)(%s)mult-row[%s]"%(in0,in1,i) for i in range(self.num_bits)]
        for rowVar in row:
            self._addBinary(rowVar, self.num_bits)
        sums = ["(%s)(%s)mult-combined[+%s]"%(in0,in1,i) for i in range(0, self.num_bits-1)] #Doesn't need to be a list of list since self.add uses full number variables, not individual bits
        sums.append("%s"%(out))

        for multiplier in range(self.num_bits): #Multiplier is the place of the bottom number (in0) for i in range(self.num_bits+multiplier): #Place within the top number (in1)
            for i in range(self.num_bits):
                if(i < multiplier): #Places
                    self._updateLinear("%s[%s]"%(row[multiplier],i), 10*scale) #Setting them to 0
                else:
                    self._and("%s[%s]"%(in0,multiplier), "%s[%s]"%(in1,i-multiplier), "%s[%s]"%(row[multiplier],i), scale=1*scale)
        for j in range(1, self.num_bits):
            if(j == 1):
                self.add(row[0], row[1], sums[1], negativeFlipAround=False, scale=1*scale)
            else:
                self.add(sums[j-1], row[j], sums[j], negativeFlipAround=False, scale=1*scale)
        for rowVar in row:
            self.reportNum(rowVar)
#        for i in range(1,len(sums)):
#            self.reportNum(sums[i])
        productBasedSignBit = "(%s)(%s)mult-productBasedSignBit"%(in0,in1) #Does not account for 0
        fullOneCheck = "(%s)(%s)mult-fullOneCheck"%(in0,in1)
        finalSignBit = "(%s)(%s)mult-finalSignBit"%(in0,in1)
        #for i in range(self.num_bits): #Creating a three-way and
        #    self._updateQuadratic("%s[%s]"%(out,i),fullOneCheck,-1*scale)
        #TODO: Create full 0 or 1 check to see if it's all 0s
        self._xor("%s[%s]"%(in0,self.NUM_BITS-1),"%s[%s]"%(in1,self.NUM_BITS-1),productBasedSignBit,scale=1*scale)
        for i in range(self.NUM_BITS, self.num_bits): #Ensuring the sign bits and all the ones before match up
            self._lockXNOR("%s[%s]"%(out,i),productBasedSignBit,scale=1*scale) #TODO: Should be finalSignBt once that's implemented
        self.num_bits = int(self.num_bits/2)

        self.reportNum(out)
    def _halfAdder(self, in0, in1, out, carry, scale=1):
        self._xor(in0, in1, out, scale=1*scale)
        self._and(in0, in1, carry, scale=1*scale)

    def _fullAdder(self, in0, in1, carryIn, out, carryOut, scale=1):
        if(not self.bqm.has_variable(in0)):
            self.addVar(in0)
        if(not self.bqm.has_variable(in1)):
            self.addVar(in1)
        if(not self.bqm.has_variable(carryIn)):
            self.addVar(carryIn)
        if(not self.bqm.has_variable(out)):
            self.addVar(out)
        if(not self.bqm.has_variable(carryOut)):
            self.addVar(carryOut)
        halfAddSum0 = "(%s)(%s)fullAdder-halfAddSum0"%(in0,in1)
        halfAddCarry0 = "(%s)(%s)fullAdder-halfAddCarry0"%(in0,in1)
        halfAddCarry1 = "(%s)(%s)fullAdder-halfAddCarry1"%(carryIn,halfAddSum0)
        self._halfAdder(in0, in1, halfAddSum0, halfAddCarry0, scale=1*scale)
        self._halfAdder(carryIn, halfAddSum0, out, halfAddCarry1, scale=1*scale)
        self._or(halfAddCarry0, halfAddCarry1, carryOut, scale=1*scale)
    def _fullAdderNoCarryOut(self, in0, in1, carryIn, out, scale=1):
        #Same thing as fullAdder but withoutt all the carryOut calculations
        if(not self.bqm.has_variable(in0)):
            self.addVar(in0)
        if(not self.bqm.has_variable(in1)):
            self.addVar(in1)
        if(not self.bqm.has_variable(carryIn)):
            self.addVar(carryIn)
        if(not self.bqm.has_variable(out)):
            self.addVar(out)
        halfAddSum0 = "(%s)(%s)fullAdder-halfAddSum0"%(in0,in1)
        self._xor(in0, in1, halfAddSum0, scale=1*scale)
        self._xor(carryIn, halfAddSum0, out, scale=1*scale)
    def _and(self, in0, in1, out, scale=1):
        if(not self.bqm.has_variable(in0)):
            self.addVar(in0)
        if(not self.bqm.has_variable(in1)):
            self.addVar(in1)
        if(not self.bqm.has_variable(out)):
            self.addVar(out)
        self._updateQuadratic(in0, in1, 1*scale)
        self._updateQuadratic(in0, out, -2*scale)
        self._updateQuadratic(in1, out, -2*scale)
        self._updateLinear(out, 3*scale)
    def _or(self, in0, in1, out, scale=1):
        if(not self.bqm.has_variable(in0)):
            self.addVar(in0)
        if(not self.bqm.has_variable(in1)):
            self.addVar(in1)
        if(not self.bqm.has_variable(out)):
            self.addVar(out)
        self._updateLinear(in0, 1*scale)
        self._updateLinear(in1, 1*scale)
        self._updateQuadratic(in0, in1, 1*scale)
        self._updateQuadratic(in0, out, -2*scale)
        self._updateQuadratic(in1, out, -2*scale)
        self._updateLinear(out, 1*scale)
    def _bothOff(self, in0, in1, out, scale=1):
        if(not self.bqm.has_variable(in0)):
            self.addVar(in0)
        if(not self.bqm.has_variable(in1)):
            self.addVar(in1)
        if(not self.bqm.has_variable(out)):
            self.addVar(out)
        auxiliary = "(%s)(%s)bothOff-or"%(in0,in1)
        self._or(in0,in1,auxiliary,scale=1*scale)
        self._not(auxiliary,out,scale=1*scale)
    def _not(self, in0, out, scale=1):
        if(not self.bqm.has_variable(in0)):
            self.addVar(in0)
        if(not self.bqm.has_variable(out)):
            self.addVar(out)
        self._updateQuadratic(in0, out, 2*scale)
        self._updateLinear(in0, -1*scale)
        self._updateLinear(out, -1*scale)
        self._updateOffset(1*scale)
    def _xnor(self, in0, in1, out, scale=1, auxScale=2): 
        if(not self.bqm.has_variable(in0)):
            self.addVar(in0)
        if(not self.bqm.has_variable(in1)):
            self.addVar(in1)
        if(not self.bqm.has_variable(out)):
            self.addVar(out)
        andAux = "(%s)(%s)and"%(in0,in1)
        self._and(in0, in1, andAux, scale=scale*auxScale)
        self._updateLinear(in0, -1*scale)
        self._updateLinear(in1, -1*scale)
        self._updateLinear(andAux, 2*scale)
        self._updateQuadratic(andAux, out, -4*scale)
        self._updateQuadratic(in0, out, 2*scale)
        self._updateQuadratic(in1, out, 2*scale)
        self._updateLinear(out, -1*scale)
    def _xor(self, in0, in1, out, scale=1, auxScale=2): 
        if(not self.bqm.has_variable(in0)):
            self.addVar(in0)
        if(not self.bqm.has_variable(in1)):
            self.addVar(in1)
        if(not self.bqm.has_variable(out)):
            self.addVar(out)
        andAux = "(%s)(%s)and"%(in0,in1)
        self._and(in0, in1, andAux, scale=scale*auxScale)
        self._updateLinear(in0, 1*scale)
        self._updateLinear(in1, 1*scale)
        self._updateLinear(andAux, -2*scale)
        self._updateQuadratic(andAux, out, 4*scale)
        self._updateQuadratic(in0, out, -2*scale)
        self._updateQuadratic(in1, out, -2*scale)
        self._updateLinear(out, 1*scale)
    def _lockXOR(self, in0, in1, scale=1):
        #Different from a regular XOR in that it doesn't have an out
        #Tries to force in0 and in1 to obey the rules of XOR, favoring in0!=in1
        if(not self.bqm.has_variable(in0)):
            self.addVar(in0)
        if(not self.bqm.has_variable(in1)):
            self.addVar(in1)
        self._updateLinear(in0, -1*scale)
        self._updateLinear(in1, -1*scale)
        self._updateQuadratic(in0, in1, 2*scale)
        self._updateOffset(1*scale)
    def _lockXNOR(self, in0, in1, scale=1):
        #Different from a regular XNOR in that it doesn't have an out
        #Tries to force in0 and in1 to obey the rules of XNOR, favoring in0=in1
        if(not self.bqm.has_variable(in0)):
            self.addVar(in0)
        if(not self.bqm.has_variable(in1)):
            self.addVar(in1)
        self._updateLinear(in0, 1*scale)
        self._updateLinear(in1, 1*scale)
        self._updateQuadratic(in0, in1, -2*scale)
    def _extend(self, num, scale=1):
        for i in range(2*self.num_bits): #Going to do a more thorough check than usual 
            if(not self.bqm.has_variable("%s[%s]"%(num, i))):
                self.addVar("%s[%s]"%(num,i))
        originalSign = "%s[%s]"%(num,self.num_bits-1)
        self._updateLinear(originalSign, self.num_bits*scale) #Same as adding 1 to scale in the following loop
        for i in range(self.num_bits, 2*self.num_bits):
            self._updateQuadratic("%s[%s]"%(num, i), originalSign, -2*scale)
            self._updateLinear("%s[%s]"%(num, i), 1*scale)
    def addNum(self, var, num, scale=10):
        #Used to add a specific number as a set of binary digits
        self.checks.append(["=", var, num])
        if(abs(num) >= 2**(self.num_bits-1)):
            raise ValueError("%s out of num_bits range"%(num))
        binary = toBinary(num, var, self.num_bits)
        for key in binary.keys():
            self.addVar(key)
            self._updateLinear(key, -scale*2*(binary[key]-0.5))
            if(binary[key] == 1):
                self._updateOffset(scale)
    def notNum(self, var, num, scale=10):
        #Used to add a specific number as a set of binary digits and ensure that variable is not that number
        self.checks.append(["!=", var, num])
        if(abs(num) >= 2**(self.num_bits-1)):
            raise ValueError("%s out of num_bits range"%(num))
        binary = toBinary(num, var, self.num_bits)
        for i, key in enumerate(binary.keys()):
            if(i > 0): #Only flipping if it's not the sign bit
                binary[key] = (1+binary[key])%2
            self.addVar(key)
            self._updateLinear(key, -scale*2*(binary[key]-0.5))
            if(binary[key] == 1):
                self._updateOffset(scale)
    def _addBinary(self, var, digits):
        #Creates unassigned binary digits
        for i in range(digits):
            self.addVar("%s[%s]"%(var,i))
    def _getNumBits(self, var):
        bits = 0
        while(self.bqm.has_variable("%s[%s]"%(var,bits))):
            bits += 1
        return bits
    def _updateLinear(self, var, bias):
        if(self.bqm.get_linear(var)):
            self.bqm.set_linear(var, self.bqm.get_linear(var)+bias)
        else:
            self.bqm.set_linear(var, bias)
    def _updateQuadratic(self, var0, var1, coupling):
        if((var0, var1) in self.bqm.quadratic.keys()):
            self.bqm.set_quadratic(var0, var1, self.bqm.get_quadratic(var0, var1)+coupling)
        else:
            self.bqm.set_quadratic(var0, var1, coupling)
    def _updateOffset(self, offset):
        #Changes the offset of the overall BQM
        #Not necessary, but useful for equations with constants
        self.bqm.add_offset(offset)
    def addVar(self, varName, linear=0):
        self.bqm.set_linear(varName, linear)
    def getBQM(self):
        return self.bqm
