import numpy as np
import matplotlib.pyplot as plt

aux = np.linspace(1, 10, 101)
dists = 1/2.5/aux

grid = np.zeros(len(dists) + 1)
for i in range(len(dists)):
    grid[i+1] = grid[i] + dists[i]

print(grid[-1])