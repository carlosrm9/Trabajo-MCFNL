import numpy as np
import matplotlib.pyplot as plt

# Constantes #

eps = 1.0
mu = 1.0
c = 1/np.sqrt(eps*mu)

# Clases #

class FDTD_Maxwell_1D:
    
    # Método de iniciacion #

    def __init__(self, dx, L=10, CFL=1.0, N=101, bounds=["PEC", "PEC"]):

        # Define el paso temporal y espacial #

        self.dx = dx
        self.dxDual = np.zeros(len(self.dx))
        self.dxDual[:-1] = (self.dx[1:] - self.dx[:-1])/2
        # self.dxDual[-1] = (self.dx[0] - self.dx[-1])/2
        self.dt = CFL * min(dx) / c

        # Inicializa los grids de acuerdo al paso espacial dado #

        self.x = np.zeros(N)
        for i in range(1,N):
            self.x[i] = self.x[i-1] + dx[i]
        self.xDual = (self.x[1:] + self.x[:-1])/2

        # Inicializa los campos eléctrico y magnético #

        self.e = np.zeros(self.x.shape)
        self.h = np.zeros(self.xDual.shape)

        # Pasa las condiciones de contorno a la clase #

        self.bounds = bounds

    # Método de paso #

    def step(self):

        # Pasa los campos eléctricos y magnéticos #

        e = self.e
        h = self.h

        # Define unas constantes útiles #

        cE = np.zeros(len(self.dxDual))
        cH = np.zeros(len(self.dx))
        cE[:] = -self.dt / self.dxDual[:] / eps
        cH[:] = self.dt / self.dx[:] / mu

        # Pasa las condiciones de contorno #

        bcL = self.bounds[0] # Izquierda
        bcR = self.bounds[1] # Derecha

        if bcL == "Mur":
            eMur = e[1]

        # Evolución del campo eléctrico #

        e[1:-1] = cE[:-1] * (h[1:] - h[:-1]) + e[1:-1]


        # Condiciones de contorno izquierdas #

        if bcL == "PEC":
            e[0] = 0.0
        elif bcL == "PMC":
            e[0] = e[0] + 2 * cE[0] * h[0]
        elif bcL == "Periodic":
            e[0] = + cE[-1] * (h[0] - h[-1]) + e[0]
        elif bcL == "Mur":
            e[0] = eMur + (c * self.dt - self.dx[0]) / (c * self.dt + self.dx[0]) * (e[1] - e[0])
        else:
            raise ValueError("Invalid Boundary Conditions on the left side")

        # Condiciones de contorno derechas #

        if bcR == "PEC":
            e[-1] = 0.0
        # elif bcR == "PMC":
        #     e[-1] = e[-1] + 2 * cE[-1] * h[-1]
        # elif bcR == "Periodic":
        #     e[-1] = + cE[-1] * (h[-1] - h[-2]) + e[-1] # Dudas en cE
        # elif bcR == "Mur":
        #     # Ni idea

        # Evolución del campo magnético #

        h[:] = cH[:] * (e[1:] - e[:-1]) + h[:]
        


    

