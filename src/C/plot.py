#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    res1 = np.genfromtxt("results_C_1.csv", delimiter=",")
    res2 = np.genfromtxt("results_C_2.csv", delimiter=",")
    
    # Do this the first time for the X component
    fig, ax = plt.subplots()
    ax.imshow(res1, cmap="viridis")
    fig.savefig("results_C_1.png")
    fig.clf()

    # Do it a second time for Y component
    fig, ax = plt.subplots()
    ax.imshow(res2, cmap="viridis")
    fig.savefig("results_C_2.png")
    fig.clf()

