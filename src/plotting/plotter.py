import statistics
import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

data = pd.read_csv("../../tests/outputs/data_2001_SED_order3_strobe4_gap2_window10.csv")
print(data.head())

def plot():
    fontdict = {'fontsize': 15}
    plt.figure(figsize=(12, 6))
    sns.regplot(x=data["edit distance"], y=data["mean"], data=data)
    plt.errorbar(x=data["edit distance"], y=data["mean"], 
        yerr=[data["mean"] - data["lower confidence bound"], data["upper confidence bound"] - data["mean"]],        fmt='o', 
        color='blue', 
        ecolor='red', 
        elinewidth=2, 
        capsize=4)
    plt.xlabel('Levenshtein Distance', fontdict=fontdict)
    plt.ylabel('Euclidean distance between strobemer frequency vectors', fontdict=fontdict)
    plt.title(f'Title (there\'s so much to say about this figure I don\'t know what to put)', fontdict=fontdict)
    plt.legend()
    plt.grid(True)
    plt.show()

plot()

