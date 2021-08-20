import csv
import random

#Designed for two alternative tasks (color and direction) and two alternative choices for each task (red or blue), (left or right)
#I_col {-1 (red), 1 (blue)}
#I_dir {-1 (left), 1 (right)}
def generate_data(n, col_trend=0.5, dir_trend=0.5, action_chance=0.5, dt=1, out="in.csv"):
    #n is the number of data points
    #Trend dictates percentage it drifts towards the correct answer (assumed 1)
    #Action chance dictates the chance of choosing a certain action. Used for populating C_col and C_dir
    with open(out, 'w') as f:
        writer = csv.writer(f)
        for i in range(n):
            I_col = 1 if random.random() < col_trend else -1
            I_dir = 1 if random.random() < dir_trend else -1
            C_col = 1 if random.random() < action_chance else 0
            C_dir = 1-C_col
            row = [I_col, I_dir, C_col, C_dir, dt]
            writer.writerow(row)

if(__name__ == "__main__"):
    generate_data(100, col_trend=0.9, dir_trend=0.8)
