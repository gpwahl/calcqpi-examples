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
    qpimap,mapheader=Read_idl('qpibulk.idl')
    layers=int(mapheader[4])
    qpifft=qpimap
    for i in range(layers):
       qpimap[i,:,:]=qpimap[i,:,:]-np.average(qpimap[i,:,:])
       qpifft[i,:,:]=np.fft.fftshift(np.abs(np.fft.fft2(qpimap[i,:,:])))
    biaslower=float(mapheader[9])
    biasupper=float(mapheader[10])
    refenergy=0.2
    nx = int(mapheader[2]) #LDOS pixels x
    ny = int(mapheader[3]) #LDOS pixels y
    xctr=nx//2
    yctr=ny//2
    vmin=0.0
    vmax=np.max(qpifft)/250.0
    
    font = {'family' : 'sans-serif',
        'sans-serif' : 'Arial',
        'size'   : 16}
    plt.rc('font', **font)

    ###############################
    #Plotting the main dispersions - horiz#
    ###############################
    fig=plt.figure(figsize=(9,6))
    gs=fig.add_gridspec(2,2,hspace=0.4,wspace=0.4)
    axs=gs.subplots()
    
    ###############################
    #Plotting the main dispersions - diag #
    ###############################
    #fig=plt.figure(figsize=(5,5))
    #axs=fig.add_subplot(111)
    #qpifftrotated=scipy.ndimage.rotate(qpifft, angle=45, axes=(2,1))
    
    #sz,sx, sy = qpifftrotated.shape
    #imgd=qpifftrotated[:,sx//2:3*sx//4,sy//2]
    imgd=qpifft[:,nx//2:3*nx//4,ny//2]
    #plt.scatter(x, EigenVals_1, c=OrbVals_1, lw=0, s=10,zorder=1)
    #axs[0].imshow(qpifft[1,:,:],extent=[0, 1.414, biaslower, biasupper],aspect='auto', cmap='gray')
    axs[0][1].imshow(imgd,extent=[0, 1, biaslower, biasupper],aspect='auto', cmap='Greys',origin='lower',vmin=vmin,vmax=vmax)
    #axs[3].imshow(qpifft[1,:,:],extent=[0, 1, biaslower, biasupper],aspect='auto', cmap='gray',vmin=0,vmax=1000)
    #axs.set_rasterized(True)
 
    axs[0][1].set_xticks([0,1])
    axs[0][1].set_yticks([-0.4,-0.2,0,0.2,0.4])
  
    #for i in range(1,6):
    #    plt.vlines(x=kpts*i, ymin=Ymin, ymax=Ymax)
    #plt.vlines(x=3 * self.kpoints, ymin=-2, ymax=2)
                #kpts=indarray[-1]
    axs[0][1].set_xlim([0.0,1.0])
    axs[0][1].set_ylim([-0.5,0.5])
    axs[0][1].text(1.0, -0.5, 'Bulk',fontsize=12, va='bottom',ha='right')
    axs[0][1].set_xlabel(r'$q_x$ ($2\pi/a$)')
    axs[0][1].set_ylabel(r'Bias ($\mathrm{V}$)')
    #plt.xticks([0, kpts,2*kpts], [r'$\overline{\Gamma}$', r'$\overline{\mathrm{M}}$', r'$\overline{\mathrm{X}}$'])
    
    axs[0][0].set_ylabel(r'$E$ ($\mathrm{eV}$)')
    axs[1][1].set_xlabel(r'$q_x$ ($2\pi/a$)')
    axs[1][0].set_xlabel(r'$k_x$ ($2\pi/a$)')
    axs[1][0].set_ylabel(r'$E$ ($\mathrm{eV}$)')
    axs[0][1].axhline(0.0, lw=0.5, color='black')

    spfmap,spfheader=Read_idl('spfbulk.idl')
    spflayers=int(spfheader[4])
   
    biaslower=float(spfheader[9])
    biasupper=float(spfheader[10])
    nx = int(spfheader[2]) #LDOS pixels x
    ny = int(spfheader[3]) #LDOS pixels y
    spfenergylayer=int((refenergy-biaslower)/(biasupper-biaslower)*spflayers)
    xctr=nx//2
    yctr=ny//2
    vmin=0.0
    vmax=np.max(spfmap)/8.0
    spfd=spfmap[:,:,yctr]
    axs[1][0].imshow(spfd,extent=[-0.5, 0.5, biaslower, biasupper],aspect='auto', cmap='Greys',origin='lower',vmin=vmin,vmax=vmax)
    #axs[3].imshow(qpifft[1,:,:],extent=[0, 1, biaslower, biasupper],aspect='auto', cmap='gray',vmin=0,vmax=1000)
    #axs.set_rasterized(True)
 
    axs[1][0].set_xticks([-0.5,0,0.5])
    #axs[1].set_yticks([-0.4,-0.2,0,0.2,0.4])
  
    #for i in range(1,6):
    #    plt.vlines(x=kpts*i, ymin=Ymin, ymax=Ymax)
    #plt.vlines(x=3 * self.kpoints, ymin=-2, ymax=2)
                #kpts=indarray[-1]
    axs[1][0].set_xlim([-0.5,0.5])
    axs[1][0].set_ylim([-0.5,0.5])
    axs[1][0].axhline(0.0, lw=0.5, color='black')
    axs[1][0].axhline(refenergy, lw=0.5, color='blue',linestyle='dashed')
    #plt.xticks([0, kpts,2*kpts], [r'$\overline{\Gamma}$', r'$\overline{\mathrm{M}}$', r'$\overline{\mathrm{X}}$'])
    
    #axs[1].set_ylabel(r'$E (\mathrm{meV})$')

    sspfmap,sspfheader=Read_idl('spfsurface.idl')
    sspflayers=int(spfheader[4])
   
    sbiaslower=float(sspfheader[9])
    sbiasupper=float(sspfheader[10])
    snx = int(sspfheader[2]) #LDOS pixels x
    sny = int(sspfheader[3]) #LDOS pixels y
    svmin=0.0
    svmax=np.max(sspfmap)/32.0

    axs[0][0].imshow(sspfmap[spfenergylayer,:,:],extent=[-0.5, 0.5, -0.5, 0.5], cmap='Oranges',origin='lower',vmin=svmin,vmax=svmax,alpha=1.0)
    axs[0][0].imshow(spfmap[spfenergylayer,:,:],extent=[-0.5, 0.5, -0.5, 0.5], cmap='Greys',origin='lower',vmin=vmin,vmax=vmax,alpha=0.8)
    axs[0][0].plot([-1,1],[0,0],c='r')
    axs[0][0].scatter([0],[0],s=2,c='k')
    axs[0][0].set_xticks([-0.5,0,0.5])
    axs[0][0].set_yticks([-0.5,0,0.5])
    axs[0][0].set_xlabel(r'$k_x$ ($2\pi/a$)')
    axs[0][0].set_ylabel(r'$k_y$ ($2\pi/a$)')
    #axs[0][0].text(0.0, 0.5, 'Surface',fontsize=12, va='top',ha='center',c='r')
    #axs[0][0].text(0.5, 0.5, 'Bulk',fontsize=12, va='top',ha='right')
    #axs[1].set_yticks([-0.4,-0.2,0,0.2,0.4])
  
    #for i in range(1,6):
    #    plt.vlines(x=kpts*i, ymin=Ymin, ymax=Ymax)
    #plt.vlines(x=3 * self.kpoints, ymin=-2, ymax=2)
                #kpts=indarray[-1]
    axs[0][0].set_xlim([-0.5,0.5])
    axs[0][0].set_ylim([-0.5,0.5])
    
    sspfd=sspfmap[:,:,yctr]
    axs[1][0].imshow(sspfd,extent=[-0.5, 0.5, sbiaslower, sbiasupper],aspect='auto', cmap='Oranges',origin='lower',vmin=svmin,vmax=svmax,alpha=0.6)
    #axs[3].imshow(qpifft[1,:,:],extent=[0, 1, biaslower, biasupper],aspect='auto', cmap='gray',vmin=0,vmax=1000)
    #axs.set_rasterized(True)
    #xs[1][0].axhline(refenergy,c='b',lw=0.5,linestyle='dashed')
    axs[1][0].set_xticks([-0.5,0,0.5])
    axs[1][0].set_yticks([-0.4,-0.2,0,0.2,0.4])
  
    #for i in range(1,6):
    #    plt.vlines(x=kpts*i, ymin=Ymin, ymax=Ymax)
    #plt.vlines(x=3 * self.kpoints, ymin=-2, ymax=2)
                #kpts=indarray[-1]
    #axs[1][0].set_xlim([-0.5,0.5])
    #axs[1][0].set_ylim([-0.5,0.5])

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
    simgd=sqpifft[:,snx//2:3*snx//4,sny//2]
    axs[1][1].imshow(simgd,extent=[0, 1, sbiaslower, sbiasupper],aspect='auto', cmap='Oranges',origin='lower',vmin=svmin,vmax=svmax)
    axs[1][1].text(1.0, -0.5, 'Surface',fontsize=12, va='bottom',ha='right')
    axs[1][1].set_xticks([0,1])
    axs[1][1].set_yticks([-0.4,-0.2,0,0.2,0.4])
  
    axs[1][1].set_xlim([0.0,1.0])
    axs[1][1].set_ylim([-0.5,0.5])
    axs[1][1].axhline(0.0, lw=0.5, color='black')
    axs[1][1].set_ylabel(r'Bias (V)')
    filename = "qpitopo.pdf"
    #plt.savefig(filename, dpi=300)
    for i in range(len(axs)):
        for j in range(len(axs[i])):
            label=f'({chr((j*len(axs)+i)+97+1)})'
            axs[i][j].text(
                0.0, 1.0, label, transform=(
                    axs[i][j].transAxes + ScaledTranslation(6/72, -24/72, fig.dpi_scale_trans)),
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
