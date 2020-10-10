# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anime

def ani_update(i):
    data = np.loadtxt('result_{:0=4}.tab'.format(i),
        skiprows=2
        )
    ax = plt.subplot(1,1,1)
    ax.cla()
    x = data[0:100,1:2]
    rho = data[0:100,2:3]
    u = data[0:100,3:4]
    p = data[0:100,4:5]
    ax.plot(x[:],rho[:,0],'r',color='red')
    ax.plot(x[:],u[:,0],'r',color='blue')
    ax.plot(x[:],p[:,0],'r',color='green')
    ax.set_ylim(0,1.25)
    ax.set_xlim(-0.5,0.5)
    ax.grid(True)
    step = str(np.round(i,decimals=1))
    ax.text(-0.1,1.35,'frame='+ step)
    plt.xlabel('x')
    plt.ylabel('rho(red), u(blue), p(green)')


#---------------------main part--------------------
fig = plt.figure()
ani = anime.FuncAnimation(fig, ani_update, interval=50, frames=4000)
ani.save('HLL_sod.mp4', writer='ffmpeg')
plt.show()
