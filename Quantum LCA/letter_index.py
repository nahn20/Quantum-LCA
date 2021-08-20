#Takes a letter and returns the corresponding base 26 alphabet index
#Uses the minimum number of characters possible
#Format is $var$_a
def letter_index(num, var):
    p = 0
    while(num >= 26**p): #Determining necessary number of digits
        p += 1
    if(num == 0):
        return "%s_a"%(var)
    out = "%s_"%(var) #Building modifiable out string
    for i in reversed(range(p)):
        digit = int(num / 26**i)
        num -= digit*26**i
        out += chr(digit+97)
    return out
