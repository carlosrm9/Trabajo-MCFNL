import pytest
import numpy as np
import matplotlib.pyplot as plt

import non_uniform as fdtd
import fdtd as fdtduniform

t0 = 0.000001
dxL = 0.1
dxR = 0.05
CFL0 = 1.0

grid = np.linspace(0, 10, 101)

grid = np.concatenate([np.arange(0, 5, 0.1), np.arange(5, 10 + 0.05, 0.05)])

def test_caja_pec_animacion():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, CFL = CFL0)

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()
        fd.animation(t=t0)

def test_caja_pmc_animacion():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, CFL = CFL0, bounds = ["PMC", "PMC"])

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()
        fd.animation(t=t0)

def test_caja_pbc_animacion():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, CFL = CFL0, bounds = ["PBC", "PBC"])

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()
        fd.animation(t=t0)

def test_mur_izq_animacion():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, CFL = CFL0, bounds = ["Mur", "PEC"])

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()
        fd.animation(t=t0)

def test_mur_der_animacion():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, CFL = CFL0, bounds = ["PEC", "Mur"])

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()
        fd.animation(t=t0)

def test_caja_pec():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, CFL = CFL0)

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()

    R = np.corrcoef(e0, fd.e)
    
    assert(R[0,1] >= 0.999)

def test_caja_pmc():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, CFL = CFL0, bounds = ["PMC", "PMC"])

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()

    R = np.corrcoef(e0, fd.e)
    
    assert(R[0,1] >= 0.999)

def test_caja_pbc():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, CFL = CFL0, bounds = ["PBC", "PBC"])

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 2.5, fd.dt):
        fd.step()

    R = np.corrcoef(np.exp(-(fd.x - 8)**2 / (2*s0**2)), fd.e)
    
    assert(R[0,1] >= 0.999)

def test_caja_mur_izq():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, CFL = CFL0, bounds = ["Mur", "PEC"])

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()

    R = np.all(fd.e < 0.00001)
    
    assert(R)

def test_caja_mur_der():
    fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, CFL = CFL0, bounds = ["PEC", "Mur"])

    x0 = 3.0; s0 = 0.75
    e0 = np.exp(-(fd.x - x0)**2 / (2*s0**2))

    fd.e[:] = e0[:]

    for i in np.arange(0, 20, fd.dt):
        fd.step()

    R = np.all(fd.e < 0.00001)
    
    assert(R)

def test_errores():

    NxRange = np.int32(np.round(np.logspace(1, 3, num=20)))
    err = np.zeros(NxRange.shape)

    for CFL in np.array([0.25, 0.5, 0.75, 1.0]):
        for i in range(len(NxRange)):

            gridLeft = np.arange(0, NxRange[i]/2, 0.1)
            gridRight = np.arange(NxRange[i]/2, NxRange[i]+0.05, 0.05)
            grid = np.concatenate([gridLeft, gridRight])

            fd = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, CFL=CFL)
            
            x0 = 3.0; s0 = 0.75
            initialField = np.exp(-(fd.x - x0)**2 / (2*s0**2))
            
            fd.e[:] = initialField[:]
            for _ in np.arange(0, 20, fd.dt):
                fd.step()
            
            finalField = fd.e    
            err[i] = np.sum(np.abs(finalField - initialField))
        plt.loglog(NxRange, err, '.-', label=CFL)
    
    plt.legend()
    plt.grid(which='both')
    plt.show()


def test_comparacion():
    gridLeft = np.arange(0, 10, dxL)
    gridRight = np.arange(10, 20+dxR, dxR)
    grid = np.concatenate([gridLeft, gridRight])

    fdnon = fdtd.FDTD_Maxwell_1D_nonuniform(x = grid, CFL = CFL0)
    fduni = fdtduniform.FDTD_Maxwell_1D(L = 20, CFL = CFL0, dx = dxL)

    for i in np.arange(0, 20, fdnon.dt):
        fdnon.step()
    for i in np.arange(0, 20, fduni.dt):
        fduni.step()
    
    plt(fdnon.x, fdnon.e)
    plt(fduni.x, fduni.e)
    plt.show