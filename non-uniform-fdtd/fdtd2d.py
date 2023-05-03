import numpy as np
import matplotlib.pyplot as plt

eps0 = 1.0
mu0 = 1.0
c0 = 1/np.sqrt(eps0*mu0)

class FDTD_Maxwell_2D():
    def __init__(self, Lx=10, Ly=10, CFL=1.0, Nx=101, Ny=101, boundaryConditions=["PEC", "PEC", "PEC", "PEC"]):
        self.x = np.linspace(0, Lx, num=Nx)
        self.xDual = (self.x[1:] + self.x[:-1])/2
        self.y = np.linspace(0, Ly, num=Ny)
        self.yDual = (self.y[1:] + self.y[:-1])/2

        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]
        self.dt = CFL * min(self.dx, self.dy) / c0

        # Modo TE
        self.Ex = np.zeros((Nx, Ny))
        self.Ey = np.zeros((Nx, Ny))
        self.Hz = np.zeros((Nx-1, Ny-1))

        self.boundaryConditions = boundaryConditions
    
    def step(self):
        Ex = self.Ex
        Ey = self.Ey
        Hz = self.Hz

        cEx = self.dt / self.dy / eps0 / 0.5
        cEy = -self.dt / self.dx / eps0 / 0.5
        cHz = -self.dt / mu0 / 0.5

        bcLx = self.boundaryConditions[0]
        bcLy = self.boundaryConditions[1]
        bcRx = self.boundaryConditions[2]
        bcRy = self.boundaryConditions[3]
       
        # Actualizamos campos eléctricos
        Ex[1:-1,1:-1] = Ex[1:-1,1:-1] + cEx * (Hz[1:,1:] + Hz[:-1,1:] - Hz[1:,:-1] - Hz[:-1,:-1])
        Ey[1:-1,1:-1] = Ey[1:-1,1:-1] + cEy * (Hz[1:,1:] + Hz[1:,:-1] - Hz[:-1,1:] - Hz[:-1,:-1])

        # Lado izquierdo
        if bcLx == "PEC":
            Ex[0,:] = 0.0
            Ey[0,:] = 0.0     
        if bcLy == "PEC":
            Ex[:,0] = 0.0
            Ey[:,0] = 0.0                          
        # elif bcLx == "PMC":
        #     Ex[0,:] = e[0,:] - 2* self.dt/self.dx/eps0*h[0,:]                  # PMC
        # elif bcL == "Periodic":
        #     e[0] =  (-self.dt / self.dx / eps0) * (h[0] - h[-1]) + e[0] # Periodica
        # elif bcL == "Mur":
        #     e[0] = eMur + (c0*self.dt-self.dx)/(c0*self.dt+self.dx)*(e[1]-e[0]) # Mur
        else:
            raise ValueError("Invalid boundary conditions on the left side")

        # Lado derecho
        if bcRx == "PEC":
            Ex[-1,:] = 0.0
            Ey[-1,:] = 0.0     
        if bcRy == "PEC":
            Ex[:,-1] = 0.0
            Ey[:,-1] = 0.0
        else:       
            raise ValueError("Invalid boundary conditions on the left side")

        Hz[:,:] = Hz[:,:] + cHz * (1/self.dx * (Ey[1:,1:] + Ey[1:,:-1] - Ey[:-1,1:] - Ey[:1,:-1]) 
                                   - 1/self.dy * (Ex[1:,1:] + Ex[:-1,1:] -Ex[1:,:-1] - Ex[:-1,:-1]))
