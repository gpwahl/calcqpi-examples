import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
from numpy import linalg as LA
import sys,math
import scipy.ndimage
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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
    qpimap,mapheader=Read_idl('qpibulk.idl')
    layers=int(mapheader[4])
    qpifft=qpimap
    for i in range(layers):
       qpimap[i,:,:]=qpimap[i,:,:]-np.average(qpimap[i,:,:])
       qpifft[i,:,:]=np.fft.fftshift(np.abs(np.fft.fft2(qpimap[i,:,:])))
    biaslower=float(mapheader[9])
    biasupper=float(mapheader[10])
    refenergy=-0.1
    nx = int(mapheader[2]) #LDOS pixels x
    ny = int(mapheader[3]) #LDOS pixels y
    xctr=nx//2
    yctr=ny//2
    vmin=0.0
    vmax=np.max(qpifft)/250.0
    
    #font = {'family' : 'sans-serif',
    #    'sans-serif' : 'Arial',
    font={'size'   : 16}
    plt.rc('font', **font)

    ###############################
    #Plotting the main dispersions - horiz#
    ###############################
    fig=plt.figure(figsize=(12,10))
    gs=fig.add_gridspec(2,3,hspace=0.4,wspace=0.4)
    axs=gs.subplots()
    axs[0, 0].axis('off')
    ###############################
    #Plotting the main dispersions - diag #
    ###############################
    #fig=plt.figure(figsize=(5,5))
    #axs=fig.add_subplot(111)
    #qpifftrotated=scipy.ndimage.rotate(qpifft, angle=45, axes=(2,1))
    
    #sz,sx, sy = qpifftrotated.shape
    #imgd=qpifftrotated[:,sx//2:3*sx//4,sy//2]
    imgd=qpifft[:,nx//2:,ny//2]
    #plt.scatter(x, EigenVals_1, c=OrbVals_1, lw=0, s=10,zorder=1)
    #axs[0].imshow(qpifft[1,:,:],extent=[0, 1.414, biaslower, biasupper],aspect='auto', cmap='gray')
    axs[0][2].imshow(imgd,extent=[0, 1, biaslower, biasupper],aspect='auto', cmap='Greys',origin='lower',vmin=vmin,vmax=vmax)
    #axs[3].imshow(qpifft[1,:,:],extent=[0, 1, biaslower, biasupper],aspect='auto', cmap='gray',vmin=0,vmax=1000)
    #axs.set_rasterized(True)
    axs[0][2].axhline(0,c='k',lw=0.5)
    axs[0][2].set_xticks([0,1])
    axs[0][2].set_yticks([-0.6,-0.3,0,0.3,0.6])
  
    #for i in range(1,6):
    #    plt.vlines(x=kpts*i, ymin=Ymin, ymax=Ymax)
    #plt.vlines(x=3 * self.kpoints, ymin=-2, ymax=2)
                #kpts=indarray[-1]
    axs[0][2].set_xlim([0.0,1.0])
    axs[0][2].set_ylim([-0.75,0.75])
    axs[0][2].set_xlabel(r'$q_x$ ($2\pi/a$)')
    axs[0][2].set_ylabel(r'Bias ($\mathrm{V}$)')
    axs[0][2].axhline(refenergy,c='b',lw=0.5,linestyle='dashed')
    axs[0][2].text(1.0, -0.75, 'Bulk',fontsize=12, va='bottom',ha='right')

    ab = inset_axes(axs[0][2], width=1.2, height=1.2, loc="upper right")
    benergylayer=int((refenergy-biaslower)/(biasupper-biaslower)*layers)
    ab.imshow(qpifft[benergylayer,:,:],extent=[-1, 1, -1, 1],cmap='Greys',origin='lower',vmin=0,vmax=0.01*np.max(qpifft[benergylayer,:,:]))
    ab.plot([0,1],[0,0],c='r')
    ab.set_xlim(-1.0,1.0)
    ab.set_ylim(-1.0,1.0)
    ab.set_xticks([-1,0,1])
    ab.set_yticks([-1,0,1])
    ab.tick_params(labelsize=10)
    #plt.xticks([0, kpts,2*kpts], [r'$\overline{\Gamma}$', r'$\overline{\mathrm{M}}$', r'$\overline{\mathrm{X}}$'])
    
    axs[0][1].set_ylabel(r'$E (\mathrm{eV})$')
    #axs[j][i].axhline(0, 0, count, lw=2, color='black')

    spfmap,spfheader=Read_idl('spfbulk.idl')
    sspfmap,sspfheader=Read_idl('spfsurface.idl')
    spflayers=int(spfheader[4])
   
    biaslower=float(spfheader[9])
    biasupper=float(spfheader[10])
    fermilayer=int((0.0-biaslower)/(biasupper-biaslower)*spflayers)
    spfenergylayer=int((refenergy-biaslower)/(biasupper-biaslower)*spflayers)
    nx = int(spfheader[2]) #LDOS pixels x
    ny = int(spfheader[3]) #LDOS pixels y
    xctr=nx//2
    yctr=ny//2
    vmin=0.0
    vmax=np.max(spfmap)/4.0

    axs[1][0].imshow(sspfmap[spfenergylayer,:,:],extent=[-0.5, 0.5, -0.5, 0.5], cmap='Oranges',origin='lower',vmin=vmin,vmax=vmax,alpha=1.0)
    axs[1][0].imshow(spfmap[spfenergylayer,:,:],extent=[-0.5, 0.5, -0.5, 0.5], cmap='Greys',origin='lower',vmin=vmin,vmax=vmax,alpha=0.8)
    axs[1][0].plot([-1,1],[0,0],c='r')
    axs[1][0].scatter([0],[0],s=2,c='k')
    axs[1][0].set_xticks([-0.5,0,0.5])
    axs[1][0].set_yticks([-0.5,0,0.5])
    axs[1][0].set_xlabel(r'$k_x$ ($2\pi/a$)')
    axs[1][0].set_ylabel(r'$k_y$ ($2\pi/a$)')
    axs[1][0].text(0.0, 0.5, 'Surface',fontsize=12, va='top',ha='center',c='r')
    axs[1][0].text(0.5, 0.5, 'Bulk',fontsize=12, va='top',ha='right')
    #axs[1].set_yticks([-0.4,-0.2,0,0.2,0.4])
  
    #for i in range(1,6):
    #    plt.vlines(x=kpts*i, ymin=Ymin, ymax=Ymax)
    #plt.vlines(x=3 * self.kpoints, ymin=-2, ymax=2)
                #kpts=indarray[-1]
    axs[1][0].set_xlim([-0.5,0.5])
    axs[1][0].set_ylim([-0.5,0.5])
    
    spfd=spfmap[:,:,yctr]
    axs[0][1].imshow(spfd,extent=[-0.5, 0.5, biaslower, biasupper],aspect='auto', cmap='Greys',origin='lower',vmin=vmin,vmax=vmax)
    axs[0][1].axhline(0,c='k',lw=0.5)
    #axs[3].imshow(qpifft[1,:,:],extent=[0, 1, biaslower, biasupper],aspect='auto', cmap='gray',vmin=0,vmax=1000)
    #axs.set_rasterized(True)
 
    axs[0][1].set_xticks([-0.5,0,0.5])
    axs[0][1].set_yticks([-0.6,-0.3,0,0.3,0.6])
    axs[0][1].axhline(refenergy,c='b',lw=0.5,linestyle='dashed')
    #for i in range(1,6):
    #    plt.vlines(x=kpts*i, ymin=Ymin, ymax=Ymax)
    #plt.vlines(x=3 * self.kpoints, ymin=-2, ymax=2)
                #kpts=indarray[-1]
    axs[0][1].set_xlim([-0.5,0.5])
    axs[0][1].set_ylim([-0.75,0.75])

    axs[0][1].set_xlabel(r'$k_x$ ($2\pi/a$)')
    axs[0][1].text(0.5, -0.75, 'Bulk',fontsize=12, va='bottom',ha='right')
    #plt.xticks([0, kpts,2*kpts], [r'$\overline{\Gamma}$', r'$\overline{\mathrm{M}}$', r'$\overline{\mathrm{X}}$'])
    
    #axs[1].set_ylabel(r'$E (\mathrm{meV})$')
  

    qpismap,smapheader=Read_idl('qpisurface.idl')
    slayers=int(smapheader[4])
    sqpifft=qpismap
    for i in range(slayers):
       qpismap[i,:,:]=qpismap[i,:,:]-np.average(qpismap[i,:,:])
       sqpifft[i,:,:]=np.fft.fftshift(np.abs(np.fft.fft2(qpismap[i,:,:])))
    sbiaslower=float(smapheader[9])
    sbiasupper=float(smapheader[10])
    snx = int(smapheader[2]) #LDOS pixels x
    sny = int(smapheader[3]) #LDOS pixels y
    sxctr=snx//2
    syctr=sny//2
    svmin=0.0
    svmax=np.max(sqpifft)/250.0
    simgd=sqpifft[:,snx//2:,sny//2]
    axs[1][2].imshow(simgd,extent=[0, 1, sbiaslower, sbiasupper],aspect='auto', cmap='Oranges',origin='lower',vmin=svmin,vmax=svmax)
 
    axs[1][2].set_xticks([0,1])
    axs[1][2].set_yticks([-0.6,-0.3,0,0.3,0.6])

    axs[1][2].set_xlim([0.0,1.0])
    axs[1][2].set_ylim([-0.75,0.75])

    axs[1][2].axhline(0,c='k',lw=0.5)

    axs[1][2].set_xlabel(r'$q_x$ ($2\pi/a$)')
    axs[1][2].set_ylabel(r'Bias ($\mathrm{V}$)')
    axs[1][2].axhline(refenergy,c='b',lw=0.5,linestyle='dashed')
    axs[1][2].text(1.0, -0.75, 'Surface',fontsize=12, va='bottom',ha='right')
    energylayer=int((refenergy-sbiaslower)/(sbiasupper-sbiaslower)*slayers)
    a = inset_axes(axs[1][2], width=1.2, height=1.2, loc="upper right")
    a.imshow(sqpifft[energylayer,:,:],extent=[-1, 1, -1, 1],cmap='Oranges',origin='lower',vmin=0,vmax=0.01*np.max(sqpifft[energylayer,:,:]))
    a.plot([0,1],[0,0],c='r')
    a.set_xlim(-1.0,1.0)
    a.set_ylim(-1.0,1.0)
    a.set_xticks([-1,0,1])
    a.set_yticks([-1,0,1])
    a.tick_params(labelsize=10)
    
    
    spfmap,spfheader=Read_idl('spfsurface.idl')
    spflayers=int(spfheader[4])
   
    biaslower=float(spfheader[9])
    biasupper=float(spfheader[10])
    nx = int(spfheader[2]) #LDOS pixels x
    ny = int(spfheader[3]) #LDOS pixels y
    xctr=nx//2
    yctr=ny//2
    vmin=0.0
    vmax=np.max(spfmap)/4.0
    spfd=spfmap[:,:,yctr]
    axs[1][1].imshow(spfd,extent=[-0.5, 0.5, biaslower, biasupper],aspect='auto', cmap='Oranges',origin='lower',vmin=vmin,vmax=vmax)
    #axs[3].imshow(qpifft[1,:,:],extent=[0, 1, biaslower, biasupper],aspect='auto', cmap='gray',vmin=0,vmax=1000)
    #axs.set_rasterized(True)
 
    axs[1][1].set_xticks([-0.5,0,0.5])
    axs[1][1].set_yticks([-0.6,-0.3,0,0.3,0.6])
  
    #for i in range(1,6):
    #    plt.vlines(x=kpts*i, ymin=Ymin, ymax=Ymax)
    #plt.vlines(x=3 * self.kpoints, ymin=-2, ymax=2)
                #kpts=indarray[-1]
    axs[1][1].set_xlim([-0.5,0.5])
    axs[1][1].set_ylim([-0.75,0.75])

    axs[1][1].axhline(0,c='k',lw=0.5)
    axs[1][1].axhline(refenergy,c='b',lw=0.5,linestyle='dashed')
    axs[1][1].set_xlabel(r'$k_x$ ($2\pi/a$)')
    axs[1][1].set_ylabel(r'$E$ ($\mathrm{eV}$)')
    axs[1][1].text(0.5, -0.75, 'Surface',fontsize=12, va='bottom',ha='right')
    filename = "qpissh.pdf"
    #plt.savefig(filename, dpi=300)
    for i in range(len(axs)):
        for j in range(len(axs[i])):
            label=f'({chr((i+j*len(axs))+97)})'
            axs[i][j].text(
                0.0, 1.0, label, transform=(
                    axs[i][j].transAxes + ScaledTranslation(4/72, -20/72, fig.dpi_scale_trans)),
                fontsize=16, va='bottom')
    plt.savefig(filename, bbox_inches='tight',transparent=True)

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
