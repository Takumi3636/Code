import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import scipy.signal

#space
xmin = 0.0
xmax = 5.8e-3
nx = 100
dx = (xmax-xmin)/nx

#time
tmin = 0.0
tmax = 4.0
nt = 20000
dt = (tmax-tmin)/nt

c = 3.0e-6 #spreading coefficient

M = 51 #maximam fourier number

Tini = 0.0 #initial temperature

INTV = 100

#initialize
F = np.zeros(nx)

T = np.linspace(tmin, tmax, nt, endpoint = False)
X = np.linspace(xmin, xmax, nx, endpoint = False)+dx/2
F = Tini
    
rsq = c*dt/(dx**2) #constant
print(rsq)

#solution
v = np.zeros((nt,nx))
Tcent = []

#initial data
v[0] = F

#figure data
fig = plt.figure()
ims = []

#simulation
for t in range(nt-1):
    ta = t*dt
    for x in range(1,nx-1):
        if x == 53:
            v[t+1,x] = (1-2*rsq)*v[t,x]+rsq*(v[t,x+1]+v[t,x-1]) + dt*1e3*0.5*(1+scipy.signal.square(0.5*np.pi*ta, duty = 0.5))
        else:
            v[t+1,x] = (1-2*rsq)*v[t,x]+rsq*(v[t,x-1]+v[t,x+1])
    #Boundary condition
    v[t+1,0] = 2*Tini - v[t+1,1]
    v[t+1,nx-1] = 2*Tini - v[t+1,nx-2]
    v[t+1,0] = Tini
    v[t+1,nx-1] = Tini
    
    counter = 0
    counter += 1
    #Updating the plt
    
    if t%INTV == 0:
                im = plt.plot(X, v[t], color="red", label="t="+str(('%.2f'%t)))
                ims.append(im)
                #print(v[t, 53])
                Tcent.append(v[t, 53])

# print(Tcent)

for x in range(1, nx-1):
      vx = v[:,x]
      vxmax = max(vx)
      tmax = np.argmax(vx)*dt
      #print(tmax)
      #print(vxmax)

newTcent = [n/max(Tcent) for n in Tcent] #normalized data
print(newTcent)

plt.title("Pulse")
plt.xlabel("x")
plt.ylabel("T(x,t)")
plt.legend()
ani = animation.ArtistAnimation(fig, ims, interval=100, repeat=False)
plt.show()
ani.save("anim.gif", writer="pillow", fps=10)

