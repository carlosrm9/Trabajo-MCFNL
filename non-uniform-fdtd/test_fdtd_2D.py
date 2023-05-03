import pytest
import numpy as np
import matplotlib.pyplot as plt
import fdtd2d

def test_pec_box():
    fd = fdtd2d.FDTD_Maxwell_2D(Lx=20,Ly=20,CFL=0.99,Nx=201,Ny=201,
                                boundaryConditions=["PEC", "PEC", "PEC", "PEC"])
            
    x0 = 3.0; y0 = 3.0 ; s0 = 0.75
    initialExfield = np.exp(-((fd.x - x0)**2 + (fd.y - y0)**2) / (2*s0**2))
    initialEyfield = np.zeros(fd.Ey.shape)
    initialHzfield = np.zeros(fd.Hz.shape)
    fd.Ex[:] = initialExfield[:]
    fd.Ey[:] = initialEyfield[:]
    fd.Hz[:] = initialHzfield[:]
    plt.figure()
    fd.step()
    plt.plot(fd.x, fd.y, fd.Ex, '.', label='Ex')
    plt.plot(fd.y, fd.y, fd.Ey, '.', label='Ex')
    plt.plot(fd.xDual, fd.yDual, fd.Hz, '.', label='Hz')
    plt.ylim(-1.1, 1.1)
    plt.xlim(fd.x[0], fd.x[-1])
    plt.grid()
    plt.figure()
    for _ in np.arange(0, 20, fd.dt):
        fd.step()
        plt.plot(fd.x, fd.e, '.')
        plt.plot(fd.xDual, fd.h, '.')
        plt.ylim(-1.1, 1.1)
        plt.xlim(fd.x[0], fd.x[-1])
        plt.grid()
        plt.pause(0.01)
        plt.cla()

