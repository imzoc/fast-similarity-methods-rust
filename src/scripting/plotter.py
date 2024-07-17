import statistics
import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns


def plot(data_file, title, x_axis_label, y_axis_label):
    data = pd.read_csv(data_file)

    fontdict = {'fontsize': 15}
    plt.figure(figsize=(12, 6))
    sns.scatterplot(x=data["edit distance"], y=data["mean"], data=data)
    plt.errorbar(x=data["edit distance"], y=data["mean"], 
        yerr=[data["mean"] - data["lower confidence bound"], data["upper confidence bound"] - data["mean"]],        fmt='o', 
        color='blue', 
        ecolor='red', 
        elinewidth=2, 
        capsize=4)
    plt.xlabel(x_axis_label, fontdict=fontdict)
    plt.ylabel(y_axis_label, fontdict=fontdict)
    plt.title(title, fontdict=fontdict)
    plt.legend()
    plt.grid(True)
    plt.show()

def strobemer():
    data_file = "../../tests/outputs/" +\
        "data_2001/" +\
        "strobemer/" +\
        "euclidean_distance/" +\
        "order_3/" +\
        "strobe_length_8/" +\
        "strobe_window_gap_5/" +\
        "strobe_window_length_10/" +\
        "step_2/" +\
        "data.csv"
    title = 'Graph (descriptive title)'
    x_axis_label = 'Levenshtein Distance'
    y_axis_label = 'Strobemer Euclidean distance'
    plot(data_file, title, x_axis_label, y_axis_label)

strobemer()