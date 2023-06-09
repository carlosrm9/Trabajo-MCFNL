import numpy as np
import matplotlib.pyplot as plt

eps0 = 1.0
mu0 = 1.0
c0 = 1/np.sqrt(eps0*mu0)

class FDTD_Maxwell_1D():
    def __init__(self, L=10, CFL=1.0, dx=0.1, boundaryConditions=["PEC", "PEC"]):
        self.x = np.arange(0, L+dx, dx)
        self.xDual = (self.x[1:] + self.x[:-1])/2

        self.dx = self.x[1] - self.x[0]
        self.dt = CFL * self.dx / c0

        self.e = np.zeros(self.x.shape)
        self.h = np.zeros(self.xDual.shape)

        self.boundaryConditions = boundaryConditions
    
    def step(self):
        e = self.e
        h = self.h

        cE = -self.dt / self.dx / eps0
        cH = -self.dt / self.dx / mu0

        bcL = self.boundaryConditions[0]
        bcR = self.boundaryConditions[1]

        if bcL == "Mur":
            eMur = e[1]

        e[1:-1] = cE * (h[1:] - h[:-1]) + e[1:-1]

        # Lado izquierdo

        if bcL == "PEC":
            e[0] = 0.0                                         # PEC
        elif bcL == "PMC":
            e[0] = e[0] - 2* self.dt/self.dx/eps0*h[0]                  # PMC
        elif bcL == "Periodic":
            e[0] =  (-self.dt / self.dx / eps0) * (h[0] - h[-1]) + e[0] # Periodica
        elif bcL == "Mur":
            e[0] = eMur + (c0*self.dt-self.dx)/(c0*self.dt+self.dx)*(e[1]-e[0]) # Mur
        else:
            raise ValueError("Invalid boundary conditions on the left side")

        # Lado derecho
        if bcL == "PEC": 
            e[-1] = 0.0

        h[:] = cH * (e[1:] - e[:-1]) + h[:]

    def animation(self, t=0.01):
        plt.plot(self.x, self.e, '*')
        plt.plot(self.xDual, self.h, '.')
        plt.ylim(-1.1, 1.1)
        plt.xlim(self.x[0], self.x[-1])
        plt.grid()
        plt.pause(t)
        plt.cla()
