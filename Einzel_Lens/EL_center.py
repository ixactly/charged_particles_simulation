# coding:utf-8

import numpy as np
import matplotlib.pyplot as plt
import pickle
from matplotlib import cm

class CartesianGrid:                                                            #基本構造
    """
        Simple class to generate a computational grid and apply boundary conditions
    """

    def __init__(self, nx=150, ny=150, nz=150, xmin=-25, xmax=25, ymin=-25, ymax=25, zmin=-50, zmax=50):
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
        
def show_potential():
    V_num = input("表示したいEinzelLens印加電圧を入力：")
    with open('Einzel_Lens\Einzel_Lens_data\einzel_lensV{}.binaryfile'.format(V_num), 'rb') as lens:
        V = pickle.load(lens)

    mesh = CartesianGrid()
    V_pre = V[:, int(mesh.ny/2), :]
    two_d_V = V_pre.transpose()

    plt.rcParams['font.family'] ='sans-serif'#使用するフォント
    plt.rcParams['xtick.direction'] = 'in'#x軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
    plt.rcParams['ytick.direction'] = 'in'#y軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
    plt.rcParams['font.size'] = 12 #フォントの大きさ
    plt.rcParams['axes.linewidth'] = 1.0# 軸の線幅edge linewidth。囲みの太さ

    x, y = np.meshgrid(mesh.z, mesh.x)
    fig, ax = plt.subplots()

    surf = ax.contourf(x, y, two_d_V, cmap=cm.coolwarm)

    plt.xlabel("x[mm]")
    plt.ylabel("y[mm]")

    plt.colorbar(surf, shrink=0.5, aspect=5).set_label("Electric Potential[V]")
    plt.gca().set_aspect('equal', adjustable='box')

    plt.show()

if __name__ == '__main__':
    show_potential()
