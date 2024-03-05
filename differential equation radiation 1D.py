import math
import numpy as np
import scipy.signal

#space
xmin = 0.0
xmax = 8.1e-3 #ワイヤー長さ(m)
nx = 100 #空間分割数
dx = (xmax-xmin)/nx #最小空間間隔
cent = 61 #ワイヤーを100分割した場合の過熱点の位置

#time
tmin = 0.0
tmax = 5.0 #時間長さ(s)
nt = 200000 #時間分割数
dt = (tmax-tmin)/nt #最小時間間隔

c = 4.0e-5 #熱拡散係数(m^2/s)

a = 5.67e-8
r = 5000

Tini = 0.0 #室温（規格化時はゼロ）
INTV = 1000 #出力するデータの時間間隔をきめるやつ

#initialize
F = np.zeros(nx) #からっぽのセル（時間ゼロのときの温度分布をいれる容器）
F = Tini #上の行の容器の中身
    
rsq = c*dt/(dx**2) #安定化条件（計算上この値が0.5以下である必要がある）
print(rsq) #安定化条件出力

#solution
v = np.zeros((nt,nx)) #からっぽのセル（全計算結果を入れる容器（行列））
Tcent = [] #からっぽのセル（過熱点の温度の時間変化を入れる用）

#initial data
v[0] = F #時間ゼロのときの温度分布指定

#simulation
for t in range(nt-1): #時間0からnt-1までnt点で以下をおこなう（vのサイズがnt*nxなので-1がついている）
    ta = t*dt #経過時間
    #Boundary condition
    v[t+1,0] = Tini #どの時間でもワイヤーの端は室温
    v[t+1,-1] = Tini #どの時間でもワイヤーの端は室温
    
    for x in range(1,nx-1): #空間0からnx-1までnx点で以下をおこなう（vのサイズがnt*nxなので-1がついている）
        if x == cent: #過熱点だけ熱拡散方程式に発熱項がつく
            v[t+1,x] = (1-2*rsq)*v[t,x]+rsq*(v[t,x+1]+v[t,x-1]) + dt*1e3*0.5*(1+scipy.signal.square(0.4*np.pi*ta, duty = 0.4)) - dt*a*r*np.pi*(xmax-xmin)*v[t,x]**4 #熱拡散方程式（後ろの項はパルス過熱項）
        else: #過熱点以外はほかの領域から流れ込んでくる熱だけ
            v[t+1,x] = (1-2*rsq)*v[t,x]+rsq*(v[t,x-1]+v[t,x+1]) - dt*a*r*np.pi*(xmax-xmin)*v[t,x]**4  #熱拡散方程式
    
    if t%INTV == 0: #出力するデータの時間間隔をきめるやつ
                Tcent.append(v[t, cent]) #過熱点のデータだけ時間順にとりだし

newTcent = [n/max(Tcent) for n in Tcent] #規格化
print(newTcent) #データ出力
