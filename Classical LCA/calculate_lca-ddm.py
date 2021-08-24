import csv
import math
from scipy.stats import norm
import matplotlib.pyplot as plt
from generate_data import generate_data
def safe_exp(x):
    #Becomes 0 if it overflows (will overflow to a tiny decimal in our case)
    #Safer to use than math.exp() and has the same effect
    try:
        return math.exp(x)
    except OverflowError:
        return 0
#Forwards process of the LCA
#Exact process might not be correct; I mostly used this to compare answers with the Quantum LCA
def calculate_drift(params, inCSV='in.csv', out='out.csv'):
    f = {'col': [], 'dir': []}
    x = {'col': [], 'dir': []}
    y = []
    W = [] #Wiener process
    rows = []
    #Uses same file format as Quantum LCA to compare answers with the same sheet
    with open('in.csv') as csv_in:
        reader = csv.reader(csv_in)
        for t,row in enumerate(reader):
            row = list(map(int, row))
            I = {'col': row[0], 'dir': row[1]}
            C = {'col': row[2], 'dir': row[3]}
            dt = row[4]
            if(t == 0):
                x_prev = {'col': 0, 'dir': 0}
                f_prev = {'col': 0, 'dir': 0}
                W_prev = 0
                y_prev = 0
            else:
                W_prev = W[t-1]
                x_prev = {'col': x['col'][t-1], 'dir': x['dir'][t-1]}
                f_prev = {'col': f['col'][t-1], 'dir': f['dir'][t-1]}
                y_prev = y[t-1]
            W.append(W_prev + norm.rvs(scale=1)) #Not sure if I'm calculating this Wiener process correctly. Based on "The Physics of Optimal Decision Making: A Formal Analysis of Models of Performance in Two-Alternative Forced-Choice Tasks". Bogacz et al. 
            x['col'].append(((1-params['lambda'])*x_prev['col']+C['col']-params['beta']*f_prev['dir'])*dt)
            x['dir'].append(((1-params['lambda'])*x_prev['dir']+C['dir']-params['beta']*f_prev['col'])*dt)
            for key in f:
                f[key].append((params['g']*x[key][t]/1)/(1+(params['g']*x[key][t]/1)**2)+0.5) #Pseudo-sigmoid used for Quantum LCA. Use to compare outputs
                #f[key].append(1/(1+safe_exp(-params['g']*x[key][t]))) #More accurate sigmoid
            drift = f['col'][t]*I['col']+f['dir'][t]*I['dir']+params['omega']*(I['col']+I['dir'])
            y.append(y_prev+drift*dt+params['w']*W[t])
            row.append(drift) #col 5
            row.append(y) #col 6
            rows.append(row)
    with open(out, 'w') as csv_out:
        writer = csv.writer(csv_out)
        for row in rows:
            writer.writerow(row)
    t = [0]
    for i in range(len(rows)-1):
        t.append(rows[i][4]+t[i])
    for action in f.keys():
        print("f %s %s: %s"%(action,i,f[action][:3]))
    for action in x.keys():
        print("x %s %s: %s"%(action,i,x[action][:3]))
    #drift = [row[5] for row in rows]
    return t, y
def graph(t, y, outImg='out.png'):
    plt.plot(t,y)
    plt.xlim(t[0], t[len(t)-1])
    maxY = 0
    for j in y:
        maxY = abs(j) if abs(j) > abs(maxY) else maxY
    plt.ylim(-maxY, maxY)
    plt.savefig(outImg)
def main():
    #generate_data(1000, col_trend=0.8, dir_trend=0.8) #Uncomment to generate new data
    params = {'omega': 0.001, 'g': 1, 'lambda': 2, 'beta': 1, 'w': 0}
    t, y = calculate_drift(params)
    graph(t, y)
if(__name__ == "__main__"):
    main()
