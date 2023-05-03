import pytest
import numpy as np
import matplotlib.pyplot as plt

import non_uniform as fdtd

t0 = 0.0001

dx0 = np.full(100, 0.1)
dx0 = np.concatenate([np.repeat(0.1, 50), np.repeat(0.01, 500)])

def test_pec_animacion():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(dx=dx0)

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()
        fd.animation(t=t0)

def test_caja_pec():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform()

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()

    R = np.corrcoef(e0, fd.e)
    
    assert(R[0,1] >= 0.999999)






