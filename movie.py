# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anime

import glob
from os.path import join, relpath, isfile

extfname  = '.tab'
path_data = "./"
filename  = "result_"

filenames = [x for x in glob.glob(join(path_data,filename+'*'+extfname))]
filenames.sort()
frames = len(filenames)

def ani_update(i):
    plt.cla()
    data  = np.loadtxt(filenames[i])
    print(filenames[i])

    x   = data[:,1:2]
    rho = data[:,2:3]
    u   = data[:,3:4]
    p   = data[:,4:5]
    
    #ax = plt.subplot(1,1,1)
    
    ax = fig.add_subplot(111)
    ax.plot(x[:],rho[:,0],'r',color='red', label='density')
    ax.plot(x[:],u[:,0],'r',color='blue', label='velocity')
    ax.plot(x[:],p[:,0],'r',color='green', label='pressure')
    ax.set_ylim(0,2.0)
    ax.set_xlim(-0.5,0.5)
    ax.grid(True)
    step = str(np.round(i,decimals=1))
    ax.text(0.0,2.1,'frame='+ step,fontsize=20)

    plt.xlabel('x')
    plt.legend(loc='upper left')


#---------------------main part--------------------
fig = plt.figure()
ani = anime.FuncAnimation(fig, ani_update, interval=50, frames=frames)
ani.save('HLL_sod.mp4', writer='ffmpeg')
plt.show()
