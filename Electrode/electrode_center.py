# coding:utf-8

import numpy as np
import matplotlib.pyplot as plt

import pickle
from matplotlib import cm

class CartesianGrid:                                                            
    """
        Simple class to generate a computational grid and apply boundary conditions
    """

    def __init__(self, nx=200, ny=200, nz=200, xmin=-21, xmax=21, ymin=-21, ymax=21, zmin=-5, zmax=25):
        self.nx, self.ny, self.nz = nx, ny, nz
        self.ntotal = nx*ny*nz

        self.xmin, self.xmax = xmin, xmax
        self.ymin, self.ymax = ymin, ymax
        self.zmin, self.zmax = zmin, zmax

        self.dx = (xmax - xmin)/(nx - 1)
        self.dy = (ymax - ymin)/(ny - 1)
        self.dz = (zmax - zmin)/(nz - 1)

        self.x = np.arange(xmin, xmax + 0.5*self.dx, self.dx)
        self.y = np.arange(ymin, ymax + 0.5*self.dy, self.dy)
        self.z = np.arange(zmin, zmax + 0.5*self.dz, self.dz)

    def create_field(self):
        return np.zeros((self.nx, self.ny), dtype=np.float)                     

    def create_meshgrid(self):                                                 
        return np.meshgrid(self.x, self.y, self.z)

def show_electrode():
    rad = int(input("引出電極の角度を入力してください："))
    a = int(input("引出電極間距離を入力してください："))
    V_init = int(input("印加電圧を入力したください：")
    with open('Electrode/Electrode_data/electrode_rad{}_a{}_V{}.binaryfile'.format(rad, a, V_init), 'rb') as electrode:
        V = pickle.load(electrode)

    mesh = CartesianGrid()
    V_pre = V[:, int(mesh.ny/2), :]
    two_d_V = V_pre.transpose()

    plt.rcParams['font.family'] ='sans-serif'
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['font.size'] = 20 
    plt.rcParams['axes.linewidth'] = 1.0

    x, y = np.meshgrid(mesh.z, mesh.y)
    fig, ax = plt.subplots()
    surf = ax.contourf(x, y, two_d_V)

    plt.xlabel("x[mm]")
    plt.ylabel("y[mm]")
    fig.colorbar(surf).set_label("Electric Potential[V]")
    plt.gca().set_aspect('equal', adjustable='box')

    plt.show()

if __name__ == "__main__":
    show_electrode()