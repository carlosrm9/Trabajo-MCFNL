import pytest
import numpy as np
import matplotlib.pyplot as plt

import non_uniform as fdtd

t0 = 0.0001

grid = np.linspace(0, 10, 101)

grid = np.concatenate([np.arange(0, 5, 0.1), np.arange(5, 10 + 0.05, 0.05)])

def test_caja_pec_animacion():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid)

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()
        fd.animation(t=t0)

def test_caja_pmc_animacion():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, bounds = ["PMC", "PMC"])

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()
        fd.animation(t=t0)

def test_caja_pbc_animacion():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, bounds = ["PBC", "PBC"])

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()
        fd.animation(t=t0)

def test_mur_izq_animacion():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, bounds = ["Mur", "PEC"])

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()
        fd.animation(t=t0)

def test_mur_der_animacion():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, bounds = ["PEC", "Mur"])

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()
        fd.animation(t=t0)

def test_caja_pec():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid)

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()

    R = np.corrcoef(e0, fd.e)
    
    assert(R[0,1] >= 0.999)

def test_caja_pmc():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, bounds = ["PMC", "PMC"])

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()

    R = np.corrcoef(e0, fd.e)
    
    assert(R[0,1] >= 0.999)

def test_caja_pbc():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, bounds = ["PBC", "PBC"])

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 2.5, fd.dt):
        fd.step()

    R = np.corrcoef(np.exp(-(fd.x - 8)**2 / (2*s0**2)), fd.e)
    
    assert(R[0,1] >= 0.999)

def test_caja_mur_izq():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, bounds = ["Mur", "PEC"])

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()

    R = np.all(fd.e < 0.00001)
    
    assert(R)

def test_caja_mur_der():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, bounds = ["PEC", "Mur"])

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()

    R = np.all(fd.e < 0.00001)
    
    assert(R)