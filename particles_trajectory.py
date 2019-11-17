# coding utf-8
from trajectory_module import trajectory1 as tra_fir
from trajectory_module import trajectory2 as tra_sec
from trajectory_module import trajectory3 as tra_thi
from trajectory_module import trajectory4 as tra_fou

import matplotlib.pyplot as plt
import numpy as np
itera = 20
sample = 21
V_einzel = int(input("アインツェルレンズ印加電圧："))
particles1 = tra_fir.calc_trajectory(itera, sample)
particles2 = tra_sec.calc_trajectory(itera, particles1[2], sample)
particles3 = tra_thi.calc_trajectory(itera, particles2[2], sample, V_einzel)
particles4 = tra_fou.calc_trajectory(itera, particles3[2], sample)

plt.rcParams['xtick.direction'] = 'in'#x軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
plt.rcParams['ytick.direction'] = 'in'#y軸の目盛線が内向き('in')か外向き('out')か双方向か('inout')
plt.rcParams['font.size'] = 12 #フォントの大きさ
plt.rcParams['axes.linewidth'] = 1.0# 軸の線幅edge linewidth。囲みの太さ

plt.figure(figsize=(11.1 ,3))

hosei = 20

theta = np.pi*float(input("電極角度θ："))/180
a = int(input("電極間距離a："))

zform1 = np.array([50, 20*np.tan(theta), 0, 2, 2+18.5*np.tan(theta), 50, 50])

yform1 = np.array([20, 20, 5, 5, 18.5, 18.5, 20])

zform2 = np.array([72, 12+13*np.tan(theta) , 12, 14, 14+11.5*np.tan(theta), 72, 72]) + a - 10
yform2 = np.array([13, 13, 5, 5, 11.5, 11.5, 13])
zform3 = np.array([82, 102, 102, 82, 82]) + hosei
zform4 = np.array([107, 137, 137, 107, 107]) + hosei
zform5 = zform4 + 35

yform3 = np.array([19.7, 19.7, 21.4, 21.4, 19.7])

zform6 = np.array([350, 385.5, 385.5, 350]) + hosei
yform4 = np.array([13, 13, -13, -13])

plt.plot(zform1, yform1, color="k")
plt.plot(zform2, yform2, color="k")
plt.plot(zform1, yform1*(-1), color="k")
plt.plot(zform2, yform2*(-1), color="k")
plt.plot(zform3, yform3, color="k")
plt.plot(zform4, yform3, color="k")
plt.plot(zform5, yform3, color="k")
plt.plot(zform3, yform3*(-1), color="k")
plt.plot(zform4, yform3*(-1), color="k")
plt.plot(zform5, yform3*(-1), color="k")
plt.plot(zform6, yform4, color="k")

for i in range(sample):
    particles1[0][i].extend(particles2[0][i])
    particles1[0][i].extend(particles3[0][i])
    particles1[0][i].extend(particles4[0][i])

    particles1[1][i].extend(particles2[1][i])
    particles1[1][i].extend(particles3[1][i])
    particles1[1][i].extend(particles4[1][i])

for i in range(int(sample/2)):
    plt.plot(particles1[0][i+1], particles1[1][i+1], color = "r")

for i in range(int((sample)/2)):
    temporary = np.array(particles1[1][i+1])
    plt.plot(particles1[0][i+1], -temporary, color = "r")

plt.xlim(-5, 415)
plt.ylim(-25, 25)
plt.xlabel("x[mm]")
plt.ylabel("y[mm]")

plt.show()