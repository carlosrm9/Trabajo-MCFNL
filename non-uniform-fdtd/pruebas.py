import numpy as np
import matplotlib.pyplot as plt

import non_uniform as fdtd


# grid = np.random.uniform(0.0, 10.0, 99)
# grid = np.sort(grid)
# grid = np.concatenate([[0.0], grid, [10.0]])

# print(grid)

# NxRange = np.int32(np.round(np.logspace(1, 3, num=20)))
# print(NxRange)
# dxRange = 10/(NxRange - 1)

# for j in range(len(NxRange)):

#     grid = np.zeros(NxRange[j])
#     for i in range(len(grid)-1):
#         grid[i+1] = (i + 0.5*np.random.uniform(-0.5, 0.5))*dxRange[j]
#     print(NxRange[j], grid[-1])
#     grid[-1] = 10

#     dx = np.zeros(len(grid)-1)
#     dx[:] = grid[1:] - grid[:-1]

#     plt.plot(grid[:-1], dx, '.-', label = NxRange[j])

# plt.legend()
# plt.show()

# print(grid[-1] - grid[-2])
# print(np.max(grid[1:] - grid[:-1]))
# print(grid)


NxRange = np.int32(np.round(np.logspace(1, 3, num=20)))
dxRange = 10/(NxRange - 1)
err = np.zeros(NxRange.shape)

for CFL in np.array([1.0]):
    for i in range(len(NxRange)):

        grid = np.zeros(NxRange[i])
        for j in range(len(grid)):
            grid[j] = (j + 0.5*np.random.uniform(-0.5, 0.5))*dxRange[i]
        grid[-1] = 10
        grid[0] = 0

        dx = np.zeros(len(grid)-1)
        dx[:] = grid[1:] - grid[:-1]

        plt.plot(grid[:-1], dx, '.-', label = NxRange[i])

        fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, CFL=CFL)
        
        #for _ in np.arange(0, 20, fd.dt):
        test=np.arange(0, 20, fd.dt)
        print(NxRange[i], grid[1] - grid[0], grid[-1] - grid[-2], fd.dt)

plt.show()
