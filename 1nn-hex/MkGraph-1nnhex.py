import sys
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import ScaledTranslation
import matplotlib.transforms as mtransforms

def Read_idl(idl_file):
    Header = []
    with open(idl_file,mode="rb") as fid:
        for i in range(12):#first 12 lines are header information
            Header.append(fid.readline())
        data= np.fromfile(fid,dtype=np.single) #rest of the binary info goes here

    nx = int(Header[2]) #LDOS pixels x
    ny = int(Header[3]) #LDOS pixels y
    layers = int(Header[4]) #number of energy layers

    data= data.reshape((layers,nx,ny)) #LDOS
    return data,Header

def do_plot(ax, Z, transform,vmin,vmax,cmap):
    arrsize=Z.shape
    nx=arrsize[0]
    ny=arrsize[1]
    im = ax.imshow(Z, interpolation='none',
                   origin='lower',
                   extent=[-nx/2, nx/2, -ny/2,ny/2], cmap=cmap,vmin=vmin, vmax=vmax,clip_on=True)

    trans_data = transform + ax.transData
    im.set_transform(trans_data)

    # display intended extent of the image
    x1, x2, y1, y2 = im.get_extent()
    ax.plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], "y--",
            transform=trans_data)
    #ax.set_xlim(-nx//2, nx//2)
    #ax.set_ylim(-ny//2, ny//2)
    ax.set_xlim(-nx/4-4, nx/4+4)
    ax.set_ylim(-ny/4-4, ny/4+4)


# prepare image and figure
#fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
Z,header = Read_idl('1nn-hex.idl')
layers = int(header[4]) #number of energy layers

qpifft=np.copy(Z)
for i in range(layers):
    qpifft[i,:,:]=Z[i,:,:]-np.average(Z[i,:,:])
    qpifft[i,:,:]=np.fft.fftshift(np.abs(np.fft.fft2(qpifft[i,:,:])))

font = {'size'   : 16}
plt.rc('font', **font)
    
# image rotation
fig=plt.figure(figsize=(6,8))
gs=fig.add_gridspec(2,2,hspace=0.4,wspace=0.3)
axs=gs.subplots()
axs[0, 0].axis('off')
layer=20
bshist,bsheader = Read_idl('spf.idl')
do_plot(axs[1][0],bshist[layer], mtransforms.Affine2D().scale(1.0, -0.866).skew_deg(30, 0),np.min(bshist[layer]),np.max(bshist[layer]),cmap='Greys')
axs[1][0].set_xlabel(r'$k_x$')
axs[1][0].set_ylabel(r'$k_y$')
axs[1][0].set_xticks([])
axs[1][0].set_yticks([])
arrsize=bshist[layer].shape
nx=arrsize[0]
ny=arrsize[1]
prefact=nx/2
xpos=np.array([0,-1.0,-1.0,0,1.0,1.0,0])
ypos=np.array([-1.154725,-0.5774,0.5774,1.154725,0.5774,-0.5774,-1.154725])
axs[1][0].set_xlim(-3.0*nx/4.0, 3.0*nx/4.0)
axs[1][0].set_ylim(-3.0*ny/4.0, 3.0*ny/4.0)
axs[1][0].plot(prefact*xpos,prefact*ypos,c='r')
layer=20
# image skew
#do_plot(ax, Z[4], mtransforms.Affine2D().scale(1.0, 0.866).skew_deg(-30, 0),0.0,0.001*np.max(Z[4]))
do_plot(axs[0][1], Z[layer], mtransforms.Affine2D().scale(1.0, 0.866).skew_deg(30, 0),np.min(Z[layer]),np.max(Z[layer]),cmap='Greys')
axs[0][1].set_xticks([])
axs[0][1].set_yticks([])
axs[0][1].set_xlabel(r'$r_x$')
axs[0][1].set_ylabel(r'$r_y$')
do_plot(axs[1][1], qpifft[layer], mtransforms.Affine2D().scale(1.0, -0.866).skew_deg(30, 0),np.min(qpifft[layer]),0.01*np.max(qpifft[layer]),cmap='Greys_r')
axs[1][1].set_xlabel(r'$q_x$')
axs[1][1].set_ylabel(r'$q_y$')
axs[1][1].set_xticks([])
axs[1][1].set_yticks([])
for i in range(len(axs)):
    for j in range(len(axs[i])):
        label=f'({chr(i*len(axs[i])+j+97)})'
        axs[i][j].text(
            0.0, 1.0, label, transform=(
                axs[i][j].transAxes + ScaledTranslation(-32/72, -8/72, fig.dpi_scale_trans)),
            fontsize=18, va='bottom')
plt.savefig('1nnhexfigure.pdf', bbox_inches='tight',transparent=True)
