import numpy as np
import matplotlib.pylab as plt

a = 100 #thermal diffusivity [mm^2/s]
length_x = 200 #mm
length_y = 100
time = 2 #seconds
nodes_x = 200
nodes_y = 100
nt = 1000
INTV = 10

#Initialization

dx = length_x / nodes_x
dy = length_y /nodes_y
dt = time / nt
Ts = 20 #initial temperature

u = np.zeros((nodes_x, nodes_y)) + Ts

rsq = a*dt/(dx*dy)

if rsq < 0.5:
    print('rsq:' + str(rsq))
else:
    print('rsq(error):' + str(rsq))

#Boundary condition

u[0, :] = 0
u[-1, :] = 100

u[:, 0] = 0
u[:, -1] = 0

# u[:, 0] = np.linspace(0, length_x, nodes_x)
# u[:, -1] = np.arange(100, 0, -100/nodes)
# u[:, -1] = np.linspace(0, length_y, nodes_y)

#Visualizing with a plot

fig, axis = plt.subplots()
pcm = axis.pcolormesh(u, cmap=plt.cm.jet, vmin=0, vmax=100)
plt.colorbar(pcm, ax=axis) #make colorvbar
axis.set_xlim([0, length_y])
axis.set_ylim([0, length_x])
axis.set_aspect('equal')

#Simulation

counter = 0

while counter < nt:
    w = u.copy() 
    for i in range(1, nodes_x-1):
        for j in range(1, nodes_y-1):
            dd_ux = (w[i+1, j] - 2*w[i, j] + w[i-1, j]) / dx**2 
            dd_uy = (w[i, j+1] - 2*w[i, j] + w[i, j-1]) / dy**2
            
            u[i, j] = dt * a * (dd_ux + dd_uy) + w[i, j]
    counter += 1
    
    #Updating the plt
    
    if counter%INTV == 0:
        pcm.set_array(u)
        axis.set_title("Distribution at t:{:.3f}[s].".format(time*counter/nt))
        plt.pause(0.01)

plt.show()
