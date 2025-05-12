import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
from numpy import linalg as LA
import sys,math
import scipy.ndimage

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

def PlotGraph():
    qpimap,mapheader=Read_idl('qpi1nnrashba.idl')
    layers=int(mapheader[4])
    qpifft=qpimap
    for i in range(layers):
       qpimap[i,:,:]=qpimap[i,:,:]-np.average(qpimap[i,:,:])
       qpifft[i,:,:]=np.fft.fftshift(np.abs(np.fft.fft2(qpimap[i,:,:])))
    biaslower=float(mapheader[9])
    biasupper=float(mapheader[10])
    showenergy=-0.2
    fermilayer=int((showenergy-biaslower)/(biasupper-biaslower)*layers)
    nx = int(mapheader[2]) #LDOS pixels x
    ny = int(mapheader[3]) #LDOS pixels y
    xctr=nx//2
    yctr=ny//2
    vmin=0.0
    vmax=0.25*np.max(qpifft[fermilayer,:,:])
    
    #'family' : 'sans-serif',
    #    'sans-serif' : 'Arial',
    font = {'size'   : 16}
    plt.rc('font', **font)

    ###############################
    #Plotting the main dispersions - horiz#
    ###############################
    fig=plt.figure(figsize=(10,8))
    gs=fig.add_gridspec(2,2,hspace=0.4,wspace=0.3)
    axs=gs.subplots()
    
    ###############################
    #Plotting the main dispersions - diag #
    ###############################
    #fig=plt.figure(figsize=(5,5))
    #axs=fig.add_subplot(111)
    #qpifftrotated=scipy.ndimage.rotate(qpifft, angle=45, axes=(2,1))

    axs[1][0].imshow(qpifft[fermilayer,:,:],extent=[-2, 2, -2, 2],cmap='Greys',origin='lower',vmin=0,vmax=vmax)

    axs[1][0].scatter([0],[0],s=2,c='k')
    axs[1][0].text(0.0, -0.08, '(0,0)',fontsize=12, va='top',ha='center')
    axs[1][0].text(1.0, -0.06, '(1,0)',fontsize=12, va='top',ha='center')
    axs[1][0].plot([-1,1],[0,0],c='r')
    #axs[1][1].axhline(0,c='k',lw=0.5)
    #axs[1][1].axvline(0,c='k',lw=0.5)
    axs[1][0].set_xticks([-1.0,0,1.0])
    axs[1][0].set_yticks([-1.0,0,1.0])
  
    axs[1][0].set_xlim([-1.5,1.5])
    axs[1][0].set_ylim([-1.5,1.5])
    axs[1][0].set_xlabel(r'$q_x$ ($2\pi/a$)')
    axs[1][0].set_ylabel(r'$q_y$ ($2\pi/a$)')
    
    #sz,sx, sy = qpifftrotated.shape
    #imgd=qpifftrotated[:,sx//2:3*sx//4,sy//2]
    imgd=qpifft[:,nx//4:3*nx//4,ny//2]
    #plt.scatter(x, EigenVals_1, c=OrbVals_1, lw=0, s=10,zorder=1)
    #axs[0].imshow(qpifft[1,:,:],extent=[0, 1.414, biaslower, biasupper],aspect='auto', cmap='gray')
    axs[1][1].imshow(imgd,extent=[-1, 1, biaslower, biasupper],aspect='auto', cmap='Greys',origin='lower',vmin=vmin,vmax=vmax)
    #axs[3].imshow(qpifft[1,:,:],extent=[0, 1, biaslower, biasupper],aspect='auto', cmap='gray',vmin=0,vmax=1000)
    #axs.set_rasterized(True)
    axs[1][1].axhline(0,c='k',lw=0.5)
    axs[1][1].axhline(showenergy,c='b',lw=0.5,linestyle='dashed')
    axs[1][1].set_xticks([-1,0,1])
    axs[1][1].set_yticks([-0.4,-0.2,0,0.2,0.4])
  
    #for i in range(1,6):
    #    plt.vlines(x=kpts*i, ymin=Ymin, ymax=Ymax)
    #plt.vlines(x=3 * self.kpoints, ymin=-2, ymax=2)
                #kpts=indarray[-1]
    axs[1][1].set_xlim([-1.0,1.0])
    axs[1][1].set_ylim([-0.5,0.1])

    #plt.xticks([0, kpts,2*kpts], [r'$\overline{\Gamma}$', r'$\overline{\mathrm{M}}$', r'$\overline{\mathrm{X}}$'])
    
    axs[1][1].set_ylabel(r'Bias ($\mathrm{V})$')
    axs[1][1].set_xlabel(r'$q_x$ ($2\pi/a$)')
    #axs[j][i].axhline(0, 0, count, lw=2, color='black')

    spfmap,spfheader=Read_idl('spf.idl')
    spflayers=int(spfheader[4])
   
    biaslower=float(spfheader[9])
    biasupper=float(spfheader[10])
    fermilayer=int((showenergy-biaslower)/(biasupper-biaslower)*layers)
    nx = int(spfheader[2]) #LDOS pixels x
    ny = int(spfheader[3]) #LDOS pixels y
    xctr=nx//2
    yctr=ny//2
    vmin=0.0
    vmax=np.max(spfmap)

    axs[0][0].imshow(spfmap[fermilayer,:,:],extent=[-0.5, 0.5, -0.5, 0.5],cmap='Greys',origin='lower',vmin=0,vmax=vmax)
    axs[0][0].scatter([0],[0],s=2,c='k')
    axs[0][0].text(0.0, -0.05, r'$\Gamma$',fontsize=12, va='top', ha='center')
    axs[0][0].text(0.0, 0.52, 'M',fontsize=12, va='bottom',ha='center')
    axs[0][0].text(0.52, 0.52, 'X',fontsize=12, va='center')
    #axs[0][0].plot([0,0],[0,0.5],c='r')
    axs[0][0].plot([-0.5,0.5],[0,0.0],c='r')
    axs[0][0].set_xticks([-0.5,0,0.5])
    axs[0][0].set_yticks([-0.5,0.0,0.5])
    axs[0][0].set_xlabel(r'$k_x$ ($2\pi/a$)')
    axs[0][0].set_ylabel(r'$k_y$ ($2\pi/a$)')
    
    spfd=spfmap[:,:,yctr]
    axs[0][1].imshow(spfd,extent=[-0.5, 0.5, biaslower, biasupper],aspect='auto', cmap='Greys',origin='lower',vmin=vmin,vmax=vmax)
    #axs[3].imshow(qpifft[1,:,:],extent=[0, 1, biaslower, biasupper],aspect='auto', cmap='gray',vmin=0,vmax=1000)
    #axs.set_rasterized(True)
 
    axs[0][1].set_xticks([-0.5,0,0.5])
    #axs[1].set_yticks([-0.4,-0.2,0,0.2,0.4])
  
    #for i in range(1,6):
    #    plt.vlines(x=kpts*i, ymin=Ymin, ymax=Ymax)
    #plt.vlines(x=3 * self.kpoints, ymin=-2, ymax=2)
                #kpts=indarray[-1]
    axs[0][1].set_xlim([-0.5,0.5])
    axs[0][1].set_ylim([-0.5,0.1])

    axs[0][1].axhline(0,c='k',lw=0.5)
    axs[0][1].axhline(showenergy,c='b',lw=0.5,linestyle='dashed')
    #plt.xticks([0, kpts,2*kpts], [r'$\overline{\Gamma}$', r'$\overline{\mathrm{M}}$', r'$\overline{\mathrm{X}}$'])
    
    axs[0][1].set_ylabel(r'$E$ ($\mathrm{eV}$)')
    axs[0][1].set_xlabel(r'$k_x$ ($2\pi/a$)')
    
    #plt.savefig(filename, dpi=300)
    for i in range(len(axs)):
        for j in range(len(axs[i])):
            label=f'({chr(i*len(axs[i])+j+97+1)})'
            axs[i][j].text(
                0.0, 1.0, label, transform=(
                    axs[i][j].transAxes + ScaledTranslation(-64/72, -8/72, fig.dpi_scale_trans)),
                fontsize=18, va='bottom')

    plt.savefig('rashba.pdf', bbox_inches='tight',transparent=True)

def showmaplayer(mapheader,qpifft,layer):
    plt.figure()
    plt.imshow(qpifft[layer,:,:],vmin=0,vmax=10)
    plt.savefig('map.pdf',transparent=True,bbox_inches='tight')

def main():
    ###To run python3 seedname Fermi_Energy Ymin Ymax
    ### e.g python3 Sr2RuO4_hr.dat 0.0 -1 1
    #cutname = str(sys.argv[1]) #Name of Hamiltonian file e.g Sr2RuO4.dat
    
    
    PlotGraph()

main()
