import pytest
import numpy as np
import matplotlib.pyplot as plt
import fdtd2d
from mpl_toolkits.mplot3d import Axes3D

def test_pec_box():
    fd = fdtd2d.FDTD_Maxwell_2D(Lx=20,Ly=20,CFL=0.99,Nx=201,Ny=201,
                                boundaryConditions=["PEC", "PEC", "PEC", "PEC"])
    X, Y = np.meshgrid(fd.x, fd.y)
    XDual, YDual = np.meshgrid(fd.xDual, fd.yDual)
    x0 = 3.0; y0 = 3.0 ; s0 = 0.75
    initialExfield = np.exp(-((X - x0)**2 + (Y - y0)**2) / (2*s0**2))
    initialEyfield = np.zeros(fd.Ey.shape)
    initialHzfield = np.zeros(fd.Hz.shape)
    fd.Ex[:] = initialExfield[:]
    fd.Ey[:] = initialEyfield[:]
    fd.Hz[:] = initialHzfield[:]
    fig = plt.figure()
    for _ in np.arange(0, 20, fd.dt):
        fd.step()
        plt.contour(XDual, YDual, fd.Hz)
        # ax.set_xlabel('x')
        # ax.set_ylabel('y')
        # ax.set_zlabel('Ex')
        plt.grid()
        plt.pause(0.01)
        plt.cla()

