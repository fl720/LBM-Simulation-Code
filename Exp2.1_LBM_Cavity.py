#coding:utf-8
import numpy as np
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.animation
import matplotlib.pyplot as plt
import time

x = 256
y = 256
interval_x = 256.0
interval_y = 256.0
rho0 = 1 
Kviscosity = 1./18
tau = 2./3
Re = 1000.0
U_top = Re * Kviscosity / interval_y 
V_top = 0 

diff_x = np.array([0, 1, 0, -1,  0, 1, -1, -1,  1])
diff_y = np.array([0, 0, 1,  0, -1, 1,  1, -1, -1])
weight = [16./36, 4./36 ,4./36 ,4./36 ,4./36 , 1./36, 1./36, 1./36, 1./36 ]

F = np.ones([x, y, 9])
F[:, :] = weight

dx = interval_x / x
dy = interval_y / y
dt = ((2 * tau - 1 ) * dx**2 ) / ( 6 * Kviscosity ) 
c = dx / dt 
XX = np.linspace(0, interval_x, x)
YY = np.linspace(0, interval_y, y)
X_grid, Y_grid = np.meshgrid(XX, YY)
X_grid = X_grid.transpose()
Y_grid = Y_grid.transpose()
u = np.zeros([x,y])
v = np.zeros([x,y])
rho = np.zeros([x,y])
F_eq = np.zeros(F.shape)

obstacle = np.ones([x, y], dtype=bool)
obstacle[ 1: -1 , 1 :] = False
obstacle[ :, -1] = False

figure = plt.figure(figsize=(1, 1), dpi = 320)
fluidImage = plt.imshow(np.ones([y, x]), origin='lower', norm=plt.Normalize(-0.05, 0.05), cmap=plt.get_cmap('bwr'), interpolation='none')
startTime = time.process_time()
frame_speed = 4

def nextFrame(arg):                            # (arg is the frame number, which we don't need)
    global startTime
    global F
    global F_eq
    global rho
    global u
    global v
    global fluidImage

    # Frame rate
    if (arg % 100 == 0) and (arg > 0):
        endTime = time.process_time()
        print ("%1.1f" % (100/(endTime-startTime)), 'frames per second')
        startTime = endTime

    for i in range(frame_speed): 
    
        # S2 streaming step
        F[ : , 1 : , 2] = F [ : , : -1 , 2] # (0, 1)
        F[ : , : -1 , 4] = F [ : , 1 : , 4] # (0, -1)
         
        F[: -1 , : , 3] = F[ 1 : , : , 3 ]# (-1, 0)
        F[ 1 :, : , 1] = F[ : -1 , : , 1 ]# (1, 0)
        
        F[1 : , : -1 , 8] = F[: -1, 1 : , 8]# (1, -1)
        F[ : -1 , 1: , 6] = F[ 1 : , : -1, 6]# (-1, 1)
        
        F[: -1 , : -1 , 7] = F[ 1 :, 1 : , 7 ]# (-1 -1)
        F[1: , 1: , 5] = F[ : -1 , : -1 ,5 ]# (1, 1)

        # Bounce back BC
        Bound = F[obstacle,:]
        Bound = Bound[:, [0,3,4,1,2,7,8,5,6]]
        F[obstacle, :] = Bound
        
        # Zou-He BC
        rho[:, -1] = ( F[ : , -1 , 0] + F[ : , -1, 1] + F[:,-1, 3] + 2 * ( F[: , -1 , 2] + F[:, -1, 5] + F[:,-1,6] ) ) / ( 1 - V_top)
        F[:,-1,4] = F[ : , -1, 2 ] + 2 * rho[:,-1] * V_top / 3
        F[: , -1 , 8 ] = F[:, -1 , 6] - 0.5 * ( F[:,-1,1] - F[:,-1,3]) - rho[:, -1] * V_top / 6 + 0.5 * rho[: , -1] * U_top  
        F[: , -1 , 7 ] = F[:, -1 , 5] + 0.5 * ( F[:,-1,1] - F[:,-1,3]) - rho[:, -1] * V_top / 6 - 0.5 * rho[: , -1] * U_top  

        # S3
        rho = np.sum(F, 2)
        
        u = np.sum(F * diff_x , 2) / rho 
        v = np.sum(F * diff_y , 2) / rho
        u[:, -1] = U_top
        v[:, -1] = V_top
        
        # S4
        for it, dx, dy, w in zip(range(9), diff_x, diff_y, weight):
            F_eq[:, :, it] = rho * w *\
                (1 + ( 3 * ( dx * u + dy * v )) / c + 9/2 * (dx * u + dy * v )**2 / c**2 - 3/2 * (u * u + v * v) / c**2 )

        # S5
        F -= (F - F_eq) / tau

    # count vorticity
    # vorticity = (np.roll(v, -1, axis=0) - np.roll(v, 1, axis=0)) / (2 * dx)\
    #     - (np.roll(u,-1,axis=1) - np.roll(u,1,axis=1)) / (2* dy)
    # vorticity = vorticity.transpose()
    # fluidImage.set_array(vorticity)
 
    # count velocity
    velocity = ( u**2 + v**2 )**0.5
    velocity = velocity.transpose()
    fluidImage.set_array(velocity)

    return fluidImage        # return the figure elements to redraw

animate = matplotlib.animation.FuncAnimation(figure, nextFrame, interval = 0)

plt.show()