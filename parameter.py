import numpy as np
import matplotlib.pyplot as plt

class InitialCondition:
    def __init__(self, H = 1e-10, m = (40e-3)/(6.02e+23), q = 1.6e-19, k = 1.38064852e-23, Te = 4, eps = 8.85418782e-12, V = 2000, d = 1, phi = 20):
        
        self.H = H          #ルンゲクッタの刻み幅
        self.m = m          #イオンの質量
        self.q = q          #電気素量
        self.k = k          #ボルツマン定数
        self.Te = Te        #電子温度
        self.eps = eps      #真空の誘電率
        self.V = V          #引出電圧
        self.d = d          #電極穴直径
        self.J = (4*eps*(V**(3/2))*np.sqrt(2*q/m))/(9*d**2)     #Chaild則による空間電荷制限電流密度
        self.phi = phi
        self.I = self.J*np.pi*d**2/4    #理論最大電流量
        self.I = 0.00005    #設定電流
        
