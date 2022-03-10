#coding:utf-8
import numpy as np
import matplotlib.pyplot as plt

x = 400
y = 100
interval_x = 0.16
interval_y = 0.04
t = 800
c = 1

Kviscosity = 1.48 * 1e-5
rho0 = 1.225

dx = interval_x / x
dy = interval_y / y
dt = dx / c

XX = np.linspace(0, interval_x, x)
YY = np.linspace(0, interval_y, y)
X_grid, Y_grid = np.meshgrid(XX, YY)
X_grid = X_grid.transpose()
Y_grid = Y_grid.transpose()

diff_x = np.array([0, 1, 0, -1,  0, 1, -1, -1,  1])
diff_y = np.array([0, 0, 1,  0, -1, 1,  1, -1, -1])
weight = [16./36, 4./36 ,4./36 ,4./36 ,4./36 , 1./36, 1./36, 1./36, 1./36 ]

F = np.ones([x, y, 9])
np.random.seed(65535)
F += 0.1 * np.random.randn(x, y, 9)
F[ : , : , 1] += 2 
rho_sum = np.sum(F,2)
for i in range(9):
    F[:, : ,i] *= rho0/rho_sum

cylinder = (X_grid - 0.04 )**2 + (Y_grid - 0.02)**2 < (0.01)**2

# tao = 3 * Kviscosity /dx + 1/2
tao = ( (6 * Kviscosity * dt) / (dx**2) + 1 ) / 2
print(tao)

rho = np.zeros([x,y])
u = np.zeros([x,y])
v = np.zeros([x,y])

fig = plt.figure(figsize=(4,1), dpi=320)

for i in range(t) :
    
    # S2
    for j in range(9):
        F[ : , : , j]  = np.roll( F[ : , : , j] , diff_y[j],  axis = 1 )
        F[ : , : , j]  = np.roll( F[ : , : , j] , diff_x[j],  axis = 0 )

    # S3
    rho = np.sum(F, 2)
    
    u = np.sum(F * diff_x , 2) / rho 
    v = np.sum(F * diff_y , 2) / rho
    
    # S4
    F_eq = np.zeros(F.shape)
    for it, dx, dy, w in zip(range(9), diff_x, diff_y, weight):
        F_eq[:, :, it] = rho * w *\
            (1 + ( 3 * ( dx * u + dy * v )) + 9/2 * (dx * u + dy * v )**2 - 3/2 * (u * u + v * v) )

    # S5
    F -= (F - F_eq) / tao

    # BC
    Bound = F[cylinder,:]
    Bound = Bound[:, [0,3,4,1,2,7,8,5,6]]
    F[cylinder, :] = Bound
    
    #Plot
    if i % 10 :
        plt.cla()
        # count vorticity
        v[cylinder] = 0
        u[cylinder] = 0
        vorticity = (np.roll(v, -1, axis=0) - np.roll(v, 1, axis=0)) / (2 * dx)\
            - (np.roll(u,-1,axis=1) - np.roll(u,1,axis=1)) / (2* dy)
        vorticity[cylinder] = 0
        vorticity = vorticity.transpose()
        plt.imshow(vorticity, cmap = 'bwr')
        plt.pause(dt)








