from amrfile import io as amrio
import numpy as np
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import os
from matplotlib.colors import LogNorm

amrio.freeAll()
#for x in range(-60000, 0, 20):
def ctof(x):
    N = len(x)+1
    xf = np.zeros(N)
    xf[1:N] = x+0.5*(x[1] - x[0])
    return xf        

def intvect(x,dx):
    return int(x[0]/dx),int(x[1]/dx)

def readplot(file, nlayer):
    
    amrID = amrio.load(file)
    time = amrio.queryTime(amrID)
    level = 0
    #first, work out the domain corners
    lo,hi = amrio.queryDomainCorners(amrID, level)
    #hi[1] = 1
    order = 0 # interpolation order, 0 for piecewise constant, 1 for linear
    x,y,thk = amrio.readBox2D(amrID, level, lo, hi, 'thickness', order)
    x,y,usrf = amrio.readBox2D(amrID, level, lo, hi, 'Z_surface', order)
    x,y,topg = amrio.readBox2D(amrID, level, lo, hi, 'Z_base', order)
    cavity = usrf - thk - topg
    x,y,ux = amrio.readBox2D(amrID, level, lo, hi, 'xVel', order)
    x,y,uy = amrio.readBox2D(amrID, level, lo, hi, 'yVel', order)

    u = np.sqrt(ux**2 + uy**2)

    sigmaf = np.array([ 0.0, 0.0712603,   0.14191442,  0.21137662,  0.27910062,  0.3445959,
                       0.40744022,  0.46728806,  0.52387456,  0.5770154,   0.62660312,  0.67260068,
                       0.71503302,  0.75397751,  0.78955403,  0.82191523,  0.85123747,  0.87771266,
                       0.90154127,  0.92292652,  0.94206966,  0.95916636,  0.97440401,  0.98795991,
                       1.])
    nlayer = len(sigmaf) - 1
    sigma = np.zeros(nlayer + 2)
    sigma[1:nlayer+1] =  0.5 * ( sigmaf[1:nlayer+1] + sigmaf[0:nlayer])
    sigma[nlayer + 1] = 1.0
    
    
    T = np.zeros( (len(y),len(x), len(sigma)))
    Tpmp = np.zeros( (len(y),len(x), len(sigma)))
    namebase = 'internalEnergy'
    for l,s in enumerate(sigma):

        name = namebase + '{:04d}'.format(l-1)
        if (l == 0):
            name = namebase + 'Surface'
        elif (l == nlayer +1 ):
            name = namebase + 'Base'
           
        #print (name,s)    
        x,y,E = amrio.readBox2D(amrID, level, lo, hi, name, order)
        E = E / 2009.0
        
        Tpmp[:,:,l] = 273.15 - 9.7456e-8 * 917.0 * 9.81 * thk * s
        T[:,:,l] = np.where(E > Tpmp[:,:,l], Tpmp[:,:,l], E)
    
    
    amrio.free(amrID)
    
    return x*1.0e-3,y*1.0e-3,time,sigma,thk,usrf,cavity,u,T,Tpmp
    

def plotframe(plot_file, frame_file):

    print (plot_file, frame_file)
    
    fig = plt.figure(figsize=(8,6))
    axL = fig.add_subplot(1,2,1)
    axR = fig.add_subplot(1,2,2,aspect='equal')

    #locations of the 'cores'
    pairs = [ (108,216+12) , (128,128), (128, 384), (192+28, 384+64)]
    colors = ['red','black','blue','pink']
    if (os.path.isfile(plot_file)):
        x,y,time,s,h,usrf,cavity,u,T,Tpmp = readplot(plot_file,24)
        theta = T - Tpmp 
        lbase = T.shape[2]-1
        dxkm = (x[1]-x[0])
        temperate_area = np.sum(np.where(theta[:,:,lbase] < -1e-10, 0, 1))
        domain_area = len(x)*len(y)
        
       
        plt.sca(axR)
        #cm = plt.pcolormesh(ctof(x),ctof(y),np.ma.masked_array(u,h<10),vmin = -1e3, vmax = 1e3, cmap = 'RdYlBu_r')
        
        cm = plt.pcolormesh(ctof(x),ctof(y),np.ma.masked_array(theta[:,:,lbase],h<10),vmin = -10, vmax = 0, cmap = 'cool')
        plt.contour(x,y,cavity,[1])
        plt.contour(x,y,h,[1])
        
        plt.contour(x,y,theta[:,:,lbase],[-1e-1])
    
        for pair,color in zip(pairs,colors):
            i,j = pair
            #TT = T[j,i,:]-273.15
            xi = x[i]
            yj = y[j]

            plt.sca(axR)
            plt.plot(xi,yj,'+',color = color)
        
            plt.sca(axL)
            z  = usrf[j,i]-s*h[j,i]
            plt.plot(T[j,i,:]-273.15,z,'.-',color=color,label='x,y = ({:3.0f},{:3.0f}) km '.format(xi,yj))
            plt.plot(Tpmp[j,i,:]-273.15,z,'--',color=color)

     
        

        plt.sca(axR)
        plt.xlabel(r'$x$ (km)')
        #plt.xlim(0,768)
        #plt.ylim(0,1382.40)


        l = len('2d.hdf5')
        step = plot_file[-(l+7):-l-1]



        plt.sca(axL)
        plt.legend()
        plt.xlabel(r'T $^\circ$C')
        plt.ylabel(r'$z$ (m)')
        #plt.xlim(-60,0)
        #plt.ylim(-500,4000)
        plt.title('time: {}'.format(int(time)))
        
                
        plt.subplots_adjust(wspace=0.15,hspace=0.15,right=0.85,bottom=0.15,left=0.15,top=0.85)


        cbar_ax = fig.add_axes([0.7, 0.95 , 0.2 , 0.025])
        cb=plt.colorbar(cm, cax=cbar_ax, 
                        orientation='horizontal',extend='min')  
        cb.set_label('(T-Tpmp)(base)')
        
        plt.savefig(frame_file,dpi=300)
        
        plt.close(fig)
        return time, temperate_area/domain_area


from glob import glob

plot_files = sorted(glob('plot*.??????.2d.hdf5'))[0::1]
N = len(plot_files)
temperate_area = np.zeros(N)
time = np.zeros(N)
for n,p in enumerate (plot_files):
    time[n],temperate_area[n] = plotframe(p,'frames/frame{:04d}.png'.format(n))

plt.figure()
plt.plot(time/1000,temperate_area,'.-')
plt.xlabel('time (ka)')
plt.ylabel('Temperate bed area fraction')
os.system('ffmpeg -framerate 12 -i frames/frame%04d.png -s:v 1280x720 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p anim.mp4')
