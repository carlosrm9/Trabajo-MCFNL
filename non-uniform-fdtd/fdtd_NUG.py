import numpy as np


eps0 = 1.0
mu0 = 1.0
c0 = 1/np.sqrt(eps0*mu0)


class FDTD_Maxwell_1D_Nonuniform_Grid():
    # para que no pasen cosas raras, Umax tiene que ser como maximo (Nx-1)/L
    def __init__(self, L = 10,CFL=1.0, Nx=101, Umin=0,Umax =0.1, boundaryConditions=["PEC", "PEC"]):
        
        site_space = np.linspace(0,L,Nx)
        self.x = site_space + np.random.uniform(Umin,Umax, size =Nx)
        self.x[0] = 0
        self.x[-1] = L
        self.xDual = (self.x[1:] + self.x[:-1])/2
        
        
        self.dx = self.x[1:] - self.x[:-1]
        self.dxDual  = (self.dx[1:] + self.dx[:-1])/2
        self.dt = CFL * min(self.dx) / c0

        self.e = np.zeros(self.x.shape)
        self.h = np.zeros(self.xDual.shape)
        
        self.boundaryConditions = boundaryConditions
        
        
    def step(self):
        
        e = self.e
        h = self.h

        cE = -self.dt / self.dxDual / eps0
        cH = self.dt / self.dx / mu0

        bcL = self.boundaryConditions[0]
        if bcL == "Mur":
            eMur = e[1]

        e[1:-1] = cE * (h[1:] - h[:-1]) + e[1:-1]

        # Lado izquierdo
        if bcL == "PEC":
            e[0] = 0.0                                
        elif bcL == "PMC":  
            e[0] = e[0] - 2* self.dt/self.dx/eps0*h[0]
        elif bcL == "Periodic":
            e[0] =  (-self.dt / self.dx / eps0) * (h[0] - h[-1]) + e[0]
        elif bcL == "Mur":
            e[0] = eMur + (c0*self.dt-self.dx)/(c0*self.dt+self.dx)*(e[1]-e[0])
        else:
            raise ValueError("Invalid boundary conditions on the left side")

        # Lado derecho
        e[-1] = 0.0
        # e[-1] = e[0]

        h[:] = cH * (e[1:] - e[:-1]) + h[:]

        print(self.x)


