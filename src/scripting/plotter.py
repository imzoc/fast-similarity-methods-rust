import statistics
import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

data = pd.read_csv("../../tests/outputs/" +
    "strobemer/" +
    "data_2001/" +
    "order_3/" + 
    "strobe_length_4/" +
    "strobe_window_gap_2/" +
    "strobe_window_length_10/" +
    "data.csv"
)
print(data.head())

def plot():
    fontdict = {'fontsize': 15}
    plt.figure(figsize=(12, 6))
    sns.scatterplot(x=data["edit distance"], y=np.log2(data["mean"]), data=data)
    plt.errorbar(x=data["edit distance"], y=np.log2(data["mean"]), 
        yerr=[data["mean"] - data["lower confidence bound"], data["upper confidence bound"] - data["mean"]],        fmt='o', 
        color='blue', 
        ecolor='red', 
        elinewidth=2, 
        capsize=4)
    plt.xlabel('Levenshtein Distance', fontdict=fontdict)
    plt.ylabel('log_2 of strobemer Euclidean distance ', fontdict=fontdict)
    plt.title(f'Title (there\'s so much to say about this figure I don\'t know what to put)', fontdict=fontdict)
    plt.legend()
    plt.grid(True)
    plt.show()

plot()

