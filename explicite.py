from math import log, exp, sin, pi
from numpy import zeros, linspace, meshgrid, ones
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def explicite(Nz, Nt, K, dz, dt):

    # Initialisation

    u = ones((Nz, Nt+1))
    profondeur = zeros(Nz+1)
    temps1D = zeros(Nt+1)

    ## Initialisation du temps de maillage

    t = linspace(0, 12, Nt)
    u = 15*ones((Nz+1, Nt+1))

    # Condition à la surface

    u[0, :] = 15-10*sin(2*pi*t/12)

    ## Boucle en temps

    maxiter = 1000
    for n in range(1, maxiter):
        uold = u
        u[:, 1] = uold[:, len(uold[0])-1]

        for i in range(1, Nt):
            for j in range(0, len(u)-1):
                profondeur[j] = u(j, i)
                if(j+1<=len(u)-1):
                    profondeur[j]+=- 2*u(j+1, i)
                if(j+2<=len(u)-1):
                    profondeur[j]+=u(j+2, i)

                profondeur[j]*=1/dz**2

                temps1D[j] = K*profondeur[j]

                for k in range(1, len(u)-2):
                    u[k, i+1] = u[k, i] + dt*temps1D[k]
                
                u[0, i+1] = u[1, i+1]
            
            if(abs(u-uold).max() < 1e-4):
                break
    return u


    ## Initialisation des paramètres numériques

Nz = 500
dz = 0.25
Nt = 5000
dt = (265*24*60*60)/Nt
K = 2*10**(-6)

u=explicite(Nz, Nt, K, dz, dt)
# plt.plot(S, C)
# plt.plot(S, C)
# plt.show()