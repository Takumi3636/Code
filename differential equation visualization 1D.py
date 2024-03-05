import numpy as np
import matplotlib.pylab as plt
import sys
import scipy.signal as signal

#property
a = 1.0e-6 #thermal diffusivity[m^2/s]

#space
length = 5.5e-3 #[m]
nodes = 100 #partitions of space
dx = length / nodes #micro section

#time
time = 2 #[s]
nt = 4000 #partitions of time
dt = time / nt

extra = 50 #pick up temp from one point
INTV = 100

rsq = a*dt/(dx**2) #stabile conditon

if rsq < 0.5:
    print('rsq:' + str(rsq))
else:
    print('rsq:' + str(rsq))
    print('error')
    sys.exit()
#initialization
Ts = 20 #initial temperature

u = np.zeros(nodes) + Ts #introduction of matrix

u[0] = 20 #boundary conditon
u[-1] = 20
Tcent = []
Tcent.append(u[extra])

#Visualization with a plot
fig, axis = plt.subplots()
pcm = axis.pcolormesh([u], cmap=plt.cm.jet, vmin=0, vmax=100)
plt.colorbar(pcm, ax=axis)
# axis.set_xlim([0, 100])
axis.set_ylim([-1, 2])
# axis.set_aspect('equal')

#simulation
counter = 0

while counter < nt:
    ta = counter*dt
    w = u.copy()
    for i in range (1, nodes-1):
        if i == extra:
            u[i] = w[i] + dt * a * (w[i-1] - 2*w[i] + w[i+1]) / dx**2 + dt*1e3*1.3*(1+signal.square(1*np.pi*ta, duty = 0.5))
        else:
            u[i] = w[i] + dt * a * (w[i-1] - 2*w[i] + w[i+1]) / dx**2
             
    counter += 1
  
# updating the plt
    if counter%INTV == 0:
        pcm.set_array([u])
        axis.set_title("Distribution at t:{:.3f}[s].".format(time*counter/nt))
        plt.pause(0.01)
        Tcent.append(u[extra])
        
        
print(Tcent)
plt.show()
