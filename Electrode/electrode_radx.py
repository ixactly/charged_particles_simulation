# coding:utf-8

import pickle
import time

import numpy as np
from scipy.sparse import csc_matrix, lil_matrix

class CartesianGrid:                                                           

    # Simple class to generate a computational grid and apply boundary conditions
    
    def __init__(self, nx=200, ny=200, nz=200, xmin=-21, xmax=21, ymin=-21, ymax=21, zmin=-5, zmax=25, rad=0, r0=5, r1=11.5, r2=18, r3=19.2):
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

        self.theta = np.pi*rad/180
        self.r0, self.r1, self.r2, self.r3 = r0, r1, r2, r3

    def create_field(self):
        return np.zeros((self.nx, self.ny), dtype=np.float)                     #条件を入れるための格納庫、配列

    def create_meshgrid(self):                                                  #軸設定、max, minで表示する範囲を決定
        return np.meshgrid(self.x, self.y, self.z)

    def set_boundary_circle(self, V_1, r_in, x_in, t=1.2):                                  #円柱状にデータを配列する                                                             #rをmeshgridに合わせて考えること
        X ,Y = np.meshgrid(self.x, self.y)
        Z = X**2 + Y**2
        r_in1 = r_in + x_in/np.tan(self.theta)
        r_in2 = r_in1 + t
        for i in range(self.nx):
            for j in range(self.ny):
                if Z[i, j] <= r_in1**2:
                    Z[i, j] = 0
                if Z[i, j] >= r_in2**2:
                    Z[i, j] = 0
                if Z[i, j] != 0:
                    Z[i, j] = V_1

        Z = Z.reshape(1, self.nx, self.ny)
        return Z

    def set_boundary_cone(self, V_1, V_2, r_in, r_out, x_in, x_out, t1, t2):
        X ,Y = np.meshgrid(self.x, self.y)
        Z = X**2 + Y**2
        r_in1 = r_in + x_in/np.tan(self.theta)               #x=0で平行になる
        r_out1 = r_out + x_out/np.tan(self.theta)
        r_in2 = r_in1 + t1
        r_out2 = r_out1 + t2

        for i in range(self.nx):
            for j in range(self.ny):
                if Z[i, j] <= r_in1**2:
                    Z[i, j] = 0
                if r_in2**2 <= Z[i, j] <= r_out1**2:
                    Z[i, j] = 0
                if Z[i, j] >= r_out2**2:
                    Z[i, j] = 0
                if r_in1**2 <= Z[i, j] <= r_in2**2:
                    Z[i, j] = V_1
                if r_out1**2 <= Z[i, j] <= r_out2**2:
                    Z[i, j] = V_2

        A = Z.reshape(1, self.nx, self.ny)

        return A
    
    def set_boundary_cone1(self, V_1, V_2, r_in, r_out, r_outer, x_in, x_out, t1, t2, t3):
        X ,Y = np.meshgrid(self.x, self.y)
        Z = X**2 + Y**2
        r_in1 = r_in + x_in/np.tan(self.theta)               #x=0で平行になる
        r_out1 = r_out + x_out/np.tan(self.theta)
        r_in2 = r_in1 + t1
        r_out2 = r_out1 + t2
        r_outer1 = r_outer
        r_outer2 = r_outer1 + t3

        for i in range(self.nx):
            for j in range(self.ny):
                if Z[i, j] <= r_in1**2:
                    Z[i, j] = 0
                if r_in2**2 <= Z[i, j] <= r_out1**2:
                    Z[i, j] = 0
                if r_out2**2 <= Z[i, j] <= r_outer1**2:
                    Z[i, j] = 0
                if Z[i, j] >= r_outer2**2:
                    Z[i, j] = 0
                if r_in1**2 <= Z[i, j] <= r_in2**2:
                    Z[i, j] = V_1
                if r_out1**2 <= Z[i, j] <= r_out2**2:
                    Z[i, j] = V_2
                if r_outer1**2 <= Z[i, j] <= r_outer2**2:
                    Z[i, j] = V_2

        A = Z.reshape(1, self.nx, self.ny)

        return A

    def make_cylinder_rad0(self, a, b, c, d, e, z1, z2, z3, z4, z5):
        Z_new = e
        for i in range(int(self.nz*(z1-self.zmin)/(self.zmax - self.zmin))):
            Z_new = np.vstack((Z_new, e))

        for i in range(int(self.nz*(z1-self.zmin)/(self.zmax - self.zmin)), int(self.nz*(z2-self.zmin)/(self.zmax - self.zmin))):
            Z_new = np.vstack((Z_new, b))

        for i in range(int(self.nz*(z2-self.zmin)/(self.zmax - self.zmin)), int(self.nz*(z3-self.zmin)/(self.zmax - self.zmin))):
            Z_new = np.vstack((Z_new, a))

        for i in range(int(self.nz*(z3-self.zmin)/(self.zmax - self.zmin)), int(self.nz*(z4-self.zmin)/(self.zmax - self.zmin))):
            Z_new = np.vstack((Z_new, b))

        for i in range(int(self.nz*(z4-self.zmin)/(self.zmax - self.zmin)), int(self.nz*(z5-self.zmin)/(self.zmax - self.zmin))):
            Z_new = np.vstack((Z_new, c))

        for i in range(int(self.nz*(z5-self.zmin)/(self.zmax - self.zmin)), int(self.nz)-1):
            Z_new = np.vstack((Z_new, d))

        return Z_new

    def set_boundary_circle_rad0(self, V_1, r1, r2):                                                                                             
        X ,Y = np.meshgrid(self.x, self.y)
        Z = X**2 + Y**2

        for i in range(self.nx):
            for j in range(self.ny):
                if Z[i, j] <= r1**2:
                    Z[i, j] = 0
                if Z[i, j] >= r2**2:
                    Z[i, j] = 0
                if Z[i, j] != 0:
                    Z[i, j] = V_1

        Z = Z.reshape(1, self.nx, self.ny)
        return Z

    def set_boundary_2circle_rad0(self, V_1, V_2, r1, r2, r3, r4):
        X ,Y = np.meshgrid(self.x, self.y)
        Z = X**2 + Y**2

        for i in range(self.nx):
            for j in range(self.ny):
                if Z[i, j] <= r1**2:
                    Z[i, j] = 0
                if r2**2 <= Z[i, j] <= r3**2:
                    Z[i, j] = 0
                if Z[i, j] >= r4**2:
                    Z[i, j] = 0
                if r1**2 <= Z[i, j] <= r2**2:
                    Z[i, j] = V_1
                if r3**2 <= Z[i, j] <= r4**2:
                    Z[i, j] = V_2
        Z = Z.reshape(1, self.nx, self.ny)
        return Z

    def make_cylinder(self, V_1, V_2, a):
        if a >= (self.r2 - self.r0)*np.tan(self.theta):
            eps = 2
            Z_new = self.set_boundary_circle(V_2, 0, 0, t=self.r3+1.2)
            t = 1.2

            for i in range(int((-1-self.zmin)/self.dz)):
                Z_new = np.vstack((Z_new, self.set_boundary_circle(V_2, 0, 0, t=self.r3+1.2)))

            for i in range(int((-1-self.zmin)/self.dz), int((0-self.zmin)/self.dz)-eps):
                Z_new = np.vstack((Z_new, self.set_boundary_circle(V_2, self.r3, 0)))

            z1 = 2
            for i in range(int((0-self.zmin)/self.dz)-eps, int((z1-self.zmin)/self.dz)-eps):
                x = (i-int((0-self.zmin)/self.dz)+eps)*self.dz
                x_out = 0
                Z_new = np.vstack((Z_new, self.set_boundary_cone(V_2, V_2, self.r0, self.r3, 0, x_out, x/np.tan(self.theta), 1.2)))

            z2 = (self.r2+t-self.r0)*np.tan(self.theta)
            for i in range(int((z1-self.zmin)/self.dz)-eps, int((z2-self.zmin)/self.dz)-eps):
                x_in = (i-int((z1-self.zmin)/self.dz)+eps)*self.dz
                x_out = 0
                Z_new = np.vstack((Z_new, self.set_boundary_cone(V_2, V_2, self.r0, self.r3, x_in, x_out, t1 = 2/np.tan(self.theta), t2=1.2)))
            
            z3 = 2 + (self.r2-self.r0)*np.tan(self.theta)
            for i in range(int((z2-self.zmin)/self.dz)-eps, int((z3-self.zmin)/self.dz)-eps):
                x_in = (i-int((z2-self.zmin)/self.dz)+eps)*self.dz
                x_out = 0
                Z_new = np.vstack((Z_new, self.set_boundary_cone(V_2, V_2, (z2-z1)/np.tan(self.theta)+self.r0, self.r3, x_in, x_out, (2-x_in)/np.tan(self.theta), 1.2)))

            z4 = z1 + a
            for i in range(int((z3-self.zmin)/self.dz)-eps, int((z4-self.zmin)/self.dz)-eps):
                x_in = 0
                x_out = 0
                Z_new = np.vstack((Z_new, self.set_boundary_circle(V_2, self.r2, x_in, t+t)))

            z5 = z4 + 2
            for i in range(int((z4-self.zmin)/self.dz)-eps, int((z5-self.zmin)/self.dz)-eps):
                x = (i-int((z4-self.zmin)/self.dz)+eps)*self.dz
                x_in = 0
                x_out = 0
                Z_new = np.vstack((Z_new, self.set_boundary_cone(V_1, V_2, self.r0, self.r2, x_in, x_out, x/np.tan(self.theta), t+1.2)))

            z6 = z4 + (self.r1-self.r0+t)*np.tan(self.theta)
            for i in range(int((z5-self.zmin)/self.dz)-eps, int((z6-self.zmin)/self.dz)-eps): 
                x_in = (i-int((z5-self.zmin)/self.dz)+eps)*self.dz
                Z_new = np.vstack((Z_new, self.set_boundary_cone(V_1, V_2, self.r0, self.r2, x_in, 0, 2/np.tan(self.theta), t+1.2)))

            z7 = z5 + (self.r1-self.r0)*np.tan(self.theta)
            for i in range(int((z6-self.zmin)/self.dz)-eps, int((z7-self.zmin)/self.dz)-eps): 
                x_in = (i-int((z6-self.zmin)/self.dz)+eps)*self.dz
                Z_new = np.vstack((Z_new, self.set_boundary_cone(V_1, V_2, (z6-z5)/np.tan(self.theta)+self.r0, self.r2, x_in, 0, (2-x_in)/np.tan(self.theta), t+1.2)))

            z8 = self.zmax
            for i in range(int((z7-self.zmin)/self.dz)-eps, int((z8-self.zmin)/self.dz)): 
                Z_new = np.vstack((Z_new, self.set_boundary_cone(V_1, V_2, self.r1, self.r2, 0, 0, t, t+1.2)))

            return Z_new

        else:
            eps = 2
            r3 = self.r2 + 2/np.tan(self.theta)
            Z_new = self.set_boundary_circle(V_2, 0, 0, t=r3+1.2)
            t = 2/np.tan(self.theta)
            for i in range(int((-1-self.zmin)/self.dz)):
                Z_new = np.vstack((Z_new, self.set_boundary_circle(V_2, 0, 0, t=r3+1.2)))

            for i in range(int((-1-self.zmin)/self.dz), int((0-self.zmin)/self.dz)-eps):
                Z_new = np.vstack((Z_new, self.set_boundary_circle(V_2, r3, 0)))

            z1 = 2
            for i in range(int((0-self.zmin)/self.dz)-eps, int((z1-self.zmin)/self.dz)-eps):
                x = (i-int((0-self.zmin)/self.dz)+eps)*self.dz
                x_out = 0
                Z_new = np.vstack((Z_new, self.set_boundary_cone(V_2, V_2, self.r0, r3, 0, x_out, x/np.tan(self.theta), 1.2)))

            z2 = z1 
            for i in range(int((z1-self.zmin)/self.dz)-eps, int((z2-self.zmin)/self.dz)-eps):
                x_in = (i-int((z1-self.zmin)/self.dz)+eps)*self.dz
                x_out = 0
                Z_new = np.vstack((Z_new, self.set_boundary_cone(V_2, V_2, self.r0, r3, x_in, x_out, t, t2=1.2)))
            
            z3 = 16
            for i in range(int((z2-self.zmin)/self.dz)-eps, int((z3-self.zmin)/self.dz)-eps):
                x_in = (i-int((z2-self.zmin)/self.dz)+eps)*self.dz
                x_out = x_in
                Z_new = np.vstack((Z_new, self.set_boundary_cone1(V_1, V_2, self.r0, (z2-z1)/np.tan(self.theta)+self.r0, r3, 0, x_out, x_in/np.tan(self.theta), t, 1.2)))

            z4 = (self.r2-self.r0)*np.tan(self.theta)+2
            for i in range(int((z3-self.zmin)/self.dz)-eps, int((z4-self.zmin)/self.dz)-eps):
                x_in = (i-int((z3-self.zmin)/self.dz)+eps)*self.dz
                x_out = x_in
                Z_new = np.vstack((Z_new, self.set_boundary_cone1(V_1, V_2, self.r0, (z3-z1)/np.tan(self.theta)+self.r0, r3, x_in, x_out, t, t, t3=1.2)))

            z5 = z3 + (self.r1-self.r0)*np.tan(self.theta)
            for i in range(int((z4-self.zmin)/self.dz)-eps, int((z5-self.zmin)/self.dz)-eps):
                x_in = (i-int((z4-self.zmin)/self.dz)+eps)*self.dz
                x_out = 0
                Z_new = np.vstack((Z_new, self.set_boundary_cone(V_1, V_2, self.r0+(z4-z3)/np.tan(self.theta), self.r2, x_in, x_out, t, t+1.2)))

            z6 = self.zmax
            for i in range(int((z5-self.zmin)/self.dz)-eps, int((z6-self.zmin)/self.dz)): 
                Z_new = np.vstack((Z_new, self.set_boundary_cone(V_1, V_2, self.r1, self.r2, 0, 0, t, t+1.2)))

            return Z_new 

    def convert_to_1d_array(self, x):
        return x.reshape(self.ntotal, 1)

    def convert_to_3d_array(self, x):
        return x.reshape(self.nz, self.ny, self.nx)

class coordinate:
    def __init__(self, i, j, k, mesh):
        self.p = i*mesh.nx*mesh.ny + j*mesh.nx + k
        self.ip1 = (i+1)*mesh.nx*mesh.ny + j*mesh.nx + k
        self.im1 = (i-1)*mesh.nx*mesh.ny + j*mesh.nx + k
        self.jp1 = i*mesh.nx*mesh.ny + (j+1)*mesh.nx + k
        self.jm1 = i*mesh.nx*mesh.ny + (j-1)*mesh.nx + k
        self.kp1 = i*mesh.nx*mesh.ny + j*mesh.nx + (k+1)
        self.km1 = i*mesh.nx*mesh.ny + j*mesh.nx + (k-1)

def calc_jacobi_matrix(mesh, Z):
    #Create sparse matrix for Jacobi method
    A = lil_matrix((mesh.ntotal, mesh.ntotal))
    i, j, k = 0, 0, 0
    for i in range(1, mesh.nz-1):
        for j in range(1, mesh.ny-1):
            for k in range(1, mesh.nx-1):
                p = coordinate(i, j, k, mesh)

                if Z[i, j, k] != 0:                                           
                    A[p.p, p.p] = 1.0

                else:
                    A[p.p, p.ip1] = 1/6
                    A[p.p, p.im1] = 1/6
                    A[p.p, p.jp1] = 1/6
                    A[p.p, p.jm1] = 1/6
                    A[p.p, p.kp1] = 1/6
                    A[p.p, p.km1] = 1/6

    return A.tocsc()

class IterationControl:

    def __init__(self, max_iter, info_interval, tolerance):
        self.max_iter = max_iter
        self.info_interval = info_interval
        self.tolerance = tolerance
        self.eps = 1.0
        self.iter = 0

    def loop(self):
        self.iter += 1
        self.output_info()

        if self.eps < self.tolerance:
            return False
        elif self.iter > self.max_iter:
            print("max iteration reached")
            return False
        else:
            return True

    def calc_epsilon(self, dx):
        self.eps = np.max(abs(dx))

    def output_info(self):
        if self.iter % self.info_interval == 0:
            print("iter = %d, eps = %.3e" % (self.iter, self.eps))

#main code
def solve_eq():

    rad = float(input("引出電極角度θ："))
    a = int(input("引出電極間距離a："))

    mesh = CartesianGrid(rad=rad)
    if rad == 0:
        A = mesh.set_boundary_circle_rad0(2000, 5, 19)
        B = mesh.set_boundary_circle_rad0(2000, 17.8, 19)
        C = mesh.set_boundary_2circle_rad0(1e-100, 2000, 5, 12.7, 17.8, 19)
        D = mesh.set_boundary_2circle_rad0(1e-100, 2000, 11.5, 12.7, 17.8, 19)
        E = mesh.set_boundary_circle_rad0(2000, 0, 19)
        electrode = mesh.make_cylinder_rad0(A, B, C, D, E, -1, 0, 2, 2+a, 4+a)

    else:
        electrode = mesh.make_cylinder(1e-10, 2000, a)

    # raplace equation

    A = calc_jacobi_matrix(mesh, electrode)

    k = mesh.convert_to_1d_array(electrode)

    iter_control = IterationControl(25000, 100, 1e-2)
    
    while iter_control.loop():
        k_new = A.dot(k)
        iter_control.calc_epsilon(k_new - k)
        k, k_new = k_new, k

    # reshape for surface plotting
    k = mesh.convert_to_3d_array(k)

    with open('Electrode/electrode_data/electrode_rad{}_a{}_remake.binaryfile'.format(rad, a), 'wb') as lens:
        pickle.dump(k, lens)

if __name__ == "__main__":
    solve_eq()
