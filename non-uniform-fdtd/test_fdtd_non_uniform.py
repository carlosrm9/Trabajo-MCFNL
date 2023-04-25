import pytest
import numpy as np
import matplotlib.pyplot as plt
import fdtd

def test_pec_box():
    fd = fdtd.FDTD_Maxwell_1D()
            
    x0 = 3.0; s0 = 0.75
    initialEfield = np.exp(-(fd.x - x0)**2 / (2*s0**2))
    initialHfield = np.exp(-(fd.xDual - x0 - fd.dx/2)**2 / (2*s0**2))
    fd.e[:] = initialEfield[:]
    fd.h[:] = initialHfield[:]
    plt.figure()
    fd.step()
    plt.plot(fd.x, fd.e, '.', label='Electric field')
    plt.plot(fd.xDual, fd.h, '.', label='Magnetic field')
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
    R = np.corrcoef(initialEfield, fd.e)
    
    assert(R[0,1] >= 0.999999)

