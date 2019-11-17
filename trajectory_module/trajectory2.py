# coding:utf-8

import numpy as np
import matplotlib.pyplot as plt
import time
import pickle

import parameter as param

class CartesianGrid:                                                            #基本構造

    def __init__(self, nx=200, ny=200, nz=200, xmin=-25, xmax=25, ymin=-25, ymax=25, zmin=25, zmax=95):
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
        return np.zeros((self.nx, self.ny), dtype=np.float)                     #条件を入れるための格納庫、配列

    def create_meshgrid(self):                                                  #軸設定、max, minで表示する範囲を決定
        return np.meshgrid(self.x, self.y, self.z)

mesh = CartesianGrid()
ic = param.InitialCondition()

H = ic.H
m = ic.m
q = ic.q
k = ic.k
Te = ic.Te
eps = ic.eps
V_extract = ic.V
d = ic.d
J = ic.J
I = ic.I

delta_y = (mesh.xmax - mesh.xmin)/(mesh.nx*1000)
delta_z = (mesh.zmax - mesh.zmin)/(mesh.nz*1000)

def Runge_Kutta(x0, a, v, h):

    k1 = v
    k2 = v+a*h/2
    k3 = v+a*h/2
    k4 = v+a*h

    x = x0 + 1000*(k1+2*k2+2*k3+k4)*h/6
    return x

def getNearestValue(zlist, ylist, value, i):                 #z方向の値を返して，対応するy方向の値に代入

    idx = np.abs(np.asarray(zlist[i]) - value).argmin()
    return ylist[i][idx]

def SpaceChargeEffect1(I, r, v, L1, L2):
    pi = np.pi
    eps1 = 8.85418782
    power = 1e12

    Q = I/v

    E = Q*power*((L1/np.sqrt(L1**2+r**2))+(L2/np.sqrt(L2**2+r**2)))/(4*pi*eps1*r*1e-3)
    return E

def SpaceChargeEffect2(I, r1, r2, v, L1, L2):
    pi = np.pi
    eps1 = 8.85418782
    power = 1e12

    Q = I/v
    E = Q*power*r1*((L1/np.sqrt(L1**2+r1**2))+(L2/np.sqrt(L2**2+r1**2)))/(4*pi*eps1*r2**2*1e-3)
    return E

def calc_trajectory(itera, data, sample):

    print("phase2")

    Ez = np.zeros((mesh.nz-1, mesh.nx-1))
    Ey = np.zeros((mesh.nz-1, mesh.nx-1))

    for a in range(itera):
        data1 = [[],[],[]]
        
        zlist = [[] for i in range(sample)]
        ylist = [[] for i in range(sample)]
        vzlist = [[] for i in range(sample)]

        Az = Ez*(q/m)
        Ay = Ey*(q/m)

        for b in range(sample):                         #軌道計算
            
            t = 0
            vz0 = data[0][b]
            vy0 = data[1][b]
            z0 = 25
            y0 = data[2][b]

            vzlist[b].append(vz0)
            zlist[b].append(z0)
            ylist[b].append(y0)

            while mesh.zmin-0.3<=z0<=mesh.zmax and -19.5<=y0<=19.5:
                
                t += H

                az = Az[int((mesh.ymax-y0)/mesh.dy), int((z0-mesh.zmin)/mesh.dz)]
                ay = Ay[int((mesh.ymax-y0)/mesh.dy), int((z0-mesh.zmin)/mesh.dz)]

                vz0 += az*H
                vy0 += ay*H

                z0 = Runge_Kutta(z0, az, vz0, H)
                y0 = Runge_Kutta(y0, ay, vy0, H)

                vzlist[b].append(vz0)
                zlist[b].append(z0)
                ylist[b].append(y0)

            data1[0].append(vz0)
            data1[1].append(vy0)
            data1[2].append(y0)

        Ez = np.zeros((mesh.nz-1, mesh.nx-1))
        Ey = np.zeros((mesh.nz-1, mesh.nx-1))

        for k in range(int(mesh.nz-1)):  #軌道から電場

            y0list = []
            vz0list = []
            for j in range(sample):
                if getNearestValue(zlist, ylist, mesh.zmin+k*mesh.dz, j) < 19.5  :
                    y0list.append(getNearestValue(zlist, ylist, mesh.zmin+k*mesh.dz, j))
                    vz0list.append(getNearestValue(zlist, vzlist, mesh.zmin+k*mesh.dz, j))

            num = len(y0list)

            r = np.nanmax(y0list)
            if r == np.inf:
                r = 4           
            v = np.mean(vz0list)

            for j in range(int((mesh.ymax-20)*mesh.ny/(mesh.ymax-mesh.ymin)), int((mesh.ymax+20)*mesh.ny/(mesh.ymax-mesh.ymin))):
                if int((mesh.ymax-r)*mesh.ny/(mesh.ymax-mesh.ymin)) <= j <= int((mesh.ymax+r)*mesh.ny/(mesh.ymax-mesh.ymin)):
                    Ey[j-1, k] = Ey[j-1, k] + SpaceChargeEffect2(I*(num/sample), (mesh.ny/2-j)*mesh.dy, r, v, k*mesh.dz+mesh.zmin, 405-(k*mesh.dz+mesh.zmin))
                    
                else:
                    Ey[j-1, k] = Ey[j-1, k] + SpaceChargeEffect1(I*(num/sample), (mesh.ny/2-j)*mesh.dy, v, k*mesh.dz+mesh.zmin, 405-(k*mesh.dz+mesh.zmin))
    print(r)    
    return zlist, ylist, data1, r