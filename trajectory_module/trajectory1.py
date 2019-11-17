# coding:utf-8

import pickle

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

import parameter as param

class CartesianGrid:                                                            #基本構造

    def __init__(self, nx=200, ny=200, nz=200, xmin=-21, xmax=21, ymin=-21, ymax=21, zmin=-5, zmax=25):
        self.nx, self.ny, self.nz = nx, ny, nz
        self.ntotal = nx*ny*nz

        self.xmin, self.xmax = xmin, xmax
        self.ymin, self.ymax = ymin, ymax
        self.zmin, self.zmax = zmin, zmax

        self.dx = (xmax - xmin)/(nx - 2)
        self.dy = (ymax - ymin)/(ny - 2)
        self.dz = (zmax - zmin)/(nz - 2)

        self.x = np.arange(xmin, xmax + 0.5*self.dx, self.dx)
        self.y = np.arange(ymin, ymax + 0.5*self.dy, self.dy)
        self.z = np.arange(zmin, zmax + 0.5*self.dz, self.dz)

    def create_field(self):
        return np.zeros((self.nx, self.ny), dtype=np.float)                     #条件を入れるための格納庫、配列

    def create_meshgrid(self):                                                  #軸設定、max, minで表示する範囲を決定
        return np.meshgrid(self.x, self.y, self.z)

rad = float(input("電極角度θ："))
a = int(input("電極間距離a："))
V = int(input("印加電圧V："))

with open('Electrode/Electrode_data/electrode_rad{}_a{}_V{}.binaryfile'.format(rad, a, V), 'rb') as lens:
    V = pickle.load(lens)

mesh = CartesianGrid()
ic = param.InitialCondition()

V_pre = V[:, int(mesh.ny/2), :]
V = V_pre.transpose()

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

z_ic = 0

Ez = np.empty((mesh.nz-1, mesh.nx-1))
Ey = np.empty((mesh.nz-1, mesh.nx-1))

delta_y = (mesh.xmax - mesh.xmin)/((mesh.nx-2)*1000)
delta_z = (mesh.zmax - mesh.zmin)/((mesh.nz-2)*1000)

for i in range(mesh.nx -1):                                 #電場リセット
    for j in range(mesh.nz -1):
        Ey[i, j] = -(V[i, j] - V[i+1, j])/delta_y
        Ez[i, j] = -(V[i, j+1] - V[i, j])/delta_z

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

def SpaceChargeEffect1(I, dx, dy, r, v, i, j):
    pi = np.pi
    eps1 = 8.85418782
    power = 1e12

    Q = I*dx*1e-3/v

    E = Q*power/(4*pi*eps1*((i*dx)**2+(j*dy)**2)*1e-6)
    return E

def SpaceChargeEffect2(I, r, v, L1, L2):
    pi = np.pi
    eps1 = 8.85418782
    power = 1e12

    Q = I/v
    E = Q*power*((L1/np.sqrt(L1**2+r**2))+(L2/np.sqrt(L2**2+r**2)))/(4*pi*eps1*r*1e-3)

    return E

def SpaceChargeEffect3(I, r1, r2, v, L1, L2):
    pi = np.pi
    eps1 = 8.85418782
    power = 1e12

    Q = I/v
    E = Q*power*r1*((L1/np.sqrt(L1*L1+r1*r1))+(L2/np.sqrt(L2*L2+r1*r1)))/(4*pi*eps1*r2*r2*1e-3)

    return E

def calc_trajectory(itera, sample):
    print("phase1")
    
    for a in range(1):
        data = [[],[],[]]
        zlist = [[] for i in range(sample)]
        ylist = [[] for i in range(sample)]
        vzlist = [[] for i in range(sample)]

        Az = Ez*(q/m)
        Ay = Ey*(q/m)
        y = np.linspace(-3.9, 3.9, sample)

        for b in range(sample):                         #軌道計算
            
            t = 0
            vz0 = 24
            vy0 = 0
            z0 = z_ic
            y0 = y[b]

            vzlist[b].append(vz0)
            zlist[b].append(z0)
            ylist[b].append(y0)

            while mesh.zmin+1<=z0<=mesh.zmax-0.001 and -11.5<=y0<=11.5:
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

            data[0].append(vz0)
            data[1].append(vy0)
            data[2].append(y0)

        Ey_tmp = np.zeros((mesh.nz-1, mesh.ny-1))
        for k in range(int((0-mesh.zmin)/mesh.dz), int((mesh.zmax-mesh.zmin)/mesh.dz)-1):  
            y0list = []
            vz0list = []
            
            for j in range(sample):
                y0list.append(getNearestValue(zlist, ylist, k*mesh.dz+mesh.zmin+z_ic, j))
                vz0list.append(getNearestValue(zlist, vzlist, k*mesh.dz+mesh.zmin+z_ic, j))
            
            r = np.max(y0list)
            v = np.mean(vz0list)

            if v <= np.sqrt(2*q*ic.phi/m):
                v = np.sqrt(2*q*ic.phi/m)

            for i in range(int((z_ic-mesh.zmin)/mesh.dz), mesh.nz-1):
                for j in range(int((-12-mesh.ymin)/mesh.dy), int(mesh.ny/2)):
                    E = SpaceChargeEffect1(I, mesh.dz, mesh.dy, r, v, i-k, (mesh.ny/2)-j)
                    Ey_tmp[j, i] = Ey_tmp[j, i] + E*(mesh.ny/2-j)*mesh.dy/np.sqrt(((i-k)*mesh.dz)**2+((mesh.ny/2-j)*mesh.dy)**2)
                    
  
    for a in range(itera):
        data = [[],[],[]]
        zlist = [[] for i in range(sample)]
        ylist = [[] for i in range(sample)]
        vzlist = [[] for i in range(sample)]

        Az = Ez*(q/m)
        Ay = Ey*(q/m)
        y = np.linspace(-4, 4, sample)

        for b in range(sample):                        
            
            t = 0
            vz0 = 0
            vy0 = 0
            z0 = z_ic
            y0 = y[b]

            vzlist[b].append(vz0)
            zlist[b].append(z0)
            ylist[b].append(y0)

            while mesh.zmin+1<=z0<=mesh.zmax-0.001 and -11.5<=y0<=11.5:
                
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

            data[0].append(vz0)
            data[1].append(vy0)
            data[2].append(y0)

        for i in range(mesh.nx -1):                                 #電場リセット
            for j in range(mesh.nz -1):
                Ey[i, j] = -(V[i, j] - V[i+1, j])/delta_y
                Ez[i, j] = -(V[i, j+1] - V[i, j])/delta_z

        for k in range(int((0-mesh.zmin)/mesh.dz), int(mesh.nz-1)):  #軌道から電場
            y0list = []
            vz0list = []
            for j in range(sample):
                y0list.append(getNearestValue(zlist, ylist, k*mesh.dz+mesh.zmin+z_ic, j))
                vz0list.append(getNearestValue(zlist, vzlist, k*mesh.dz+mesh.zmin+z_ic, j))
            
            r = np.nanmax(y0list)
            if r == np.inf:
                r = 4
            v = np.mean(vz0list)
            r1 = 4

            Er = Ey_tmp[int((mesh.ymax-r)/mesh.dy)-1, k]

            for j in range(int((mesh.ymax-12)/mesh.dy), int(mesh.ny/2)):
                
                if int((mesh.ymax-r1)/mesh.dy) <= j <= int((mesh.ymax+r1)/mesh.dy): 
                    Ey[j-1, k] = Ey[j-1, k] + SpaceChargeEffect2(I, (mesh.ny/2-j)*mesh.dy, np.sqrt(2*q*V_extract/m), k*mesh.dz-mesh.zmax, 405-(k*mesh.dz-mesh.zmax))*((mesh.ny/2-j)*mesh.dy)/r1 
                else:
                    Ey[j-1, k] = Ey[j-1, k] + SpaceChargeEffect2(I, (mesh.ny/2-j)*mesh.dy, np.sqrt(2*q*V_extract/m), k*mesh.dz-mesh.zmax, 405-(k*mesh.dz-mesh.zmax)) 

                if int((mesh.ymax-r)/mesh.dy) <= j <= int((mesh.ymax+r)/mesh.dy):
                    Ey[j-1, k] = Ey[j-1, k] + Er*((mesh.ny/2-j)*mesh.dy/r)
                else:
                    Ey[j-1, k] = Ey[j-1, k] + Er*r/(((mesh.ny/2-j)*mesh.dy))

            for j in range(int(mesh.ny/2)+1, int((mesh.ymax+12)/mesh.dy)):
                Er = Ey_tmp[int((mesh.ymax-r)/mesh.dy)-1, k]
                
                if int((mesh.ymax-r1)/mesh.dy) <= j <= int((mesh.ymax+r1)/mesh.dy): 
                    Ey[j-1, k] = Ey[j-1, k] + SpaceChargeEffect2(I, (mesh.ny/2-j)*mesh.dy, np.sqrt(2*q*V_extract/m), k*mesh.dz-mesh.zmax, 405-(k*mesh.dz-mesh.zmax))*((mesh.ny/2-j)*mesh.dy)/r1 
                else:
                    Ey[j-1, k] = Ey[j-1, k] + SpaceChargeEffect2(I, (mesh.ny/2-j)*mesh.dy, np.sqrt(2*q*V_extract/m), k*mesh.dz-mesh.zmax, 405-(k*mesh.dz-mesh.zmax)) 

                if int((mesh.ymax-r)/mesh.dy) <= j <= int((mesh.ymax+r)/mesh.dy):
                    Ey[j-1, k] = Ey[j-1, k] + Er*((mesh.ny/2-j)*mesh.dy/r)
                else:
                    Ey[j-1, k] = Ey[j-1, k] + Er*r/(((mesh.ny/2-j)*mesh.dy))

    return zlist, ylist, data