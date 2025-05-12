import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
from numpy import linalg as LA
import sys,math
import scipy.ndimage
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

def loadcut(cutname):
    cutarray=np.loadtxt(cutname)
    count = len(cutarray)
    indarray=cutarray[:,0:1]
    indices=indarray[-1]
    kvectors=cutarray[:,1:4]
    energies=cutarray[:,4:5]
    orbvals=cutarray[:,5:]
    energies=energies
    store=[]
    spinstore=[]
    for x in range(len(energies)):
        Orb_dxz = 0
        Orb_dyz = 0
        Orb_dxy = 0
        spin=0
        for i in range(int(len(orbvals[0])/3)):
            Orb_dxz += abs(orbvals[x][i*3 +0])
            Orb_dyz += abs(orbvals[x][i*3 +1])
            Orb_dxy += abs(orbvals[x][i*3+ 2])
        colour=(int(Orb_dxz*0xff0000)&0xff0000)+(int(Orb_dyz*0xff00)&0xff00)+(int(Orb_dxy*0xff)&0xff)
        #print("#%06X"%colour)
        store.append("#%06X"%colour)
        spin=0
        spinno=int(len(orbvals[0])/2)
        for i in range(spinno):
            spin += abs(orbvals[x][i])
            spin -= abs(orbvals[x][spinno+i])
        spinstore.append(spin)
        
        #print(Orb_dxy, Orb_dxz, Orb_dyz)
        

    indarray=indarray/indices

    return indarray,energies,store,spinstore

def MkPlot(uspmapname,qpifftmapname):
    uspmap,uspmapheader=Read_idl(uspmapname)
 
    uspbiaslower=float(uspmapheader[9])
    uspbiasupper=float(uspmapheader[10])
    
    
    #font = {'family' : 'sans-serif',
    #    'sans-serif' : 'Arial',
    font={'size'   : 16}
    plt.rc('font', **font)

    ###############################
    #Plotting the main dispersions - horiz#
    ###############################
    fig=plt.figure(figsize=(8,5))
    gs=fig.add_gridspec(2,2)
    #gs=fig.add_gridspec(2,3, hspace=0,wspace=0)

    axs=gs.subplots()
    #(sharex='col', sharey='row')
  
    uspnx = int(uspmapheader[2]) #LDOS pixels x
    uspny = int(uspmapheader[3]) #LDOS pixels y
    layers=int(uspmapheader[4])
    

    fermilayer=int((0.0-uspbiaslower)/(uspbiasupper-uspbiaslower)*layers)
    
    axs[0][0].imshow(uspmap[fermilayer,:,:],extent=[-0.5, 0.5, -0.5, 0.5],cmap='Greys',origin='lower',vmin=0,vmax=np.max(uspmap[fermilayer,:,:]))
    axs[0][0].scatter([0],[0],s=2,c='k')
    axs[0][0].text(0.0, -0.05, r'$\Gamma$',fontsize=12, va='top', ha='center')
    axs[0][0].text(0.0, 0.52, 'M',fontsize=12, va='bottom',ha='center')
    axs[0][0].text(0.52, 0.52, 'X',fontsize=12, va='center')
    axs[0][0].plot([0,0],[0,0.5],c='r')
    axs[0][0].plot([0,0.5],[0,0.5],c='r')
    axs[0][0].set_xticks([-0.5,0,0.5])
    axs[0][0].set_yticks([-0.5,0.0,0.5])
    axs[0][0].set_xlabel(r'$k_x$ ($2\pi/a$)')
    axs[0][0].set_ylabel(r'$k_y$ ($2\pi/a$)')
    #axs[0][0].axhline(0,c='k',lw=0.5)
    #axs[0][0].axvline(0,c='k',lw=0.5)
    
    imgh=uspmap[:,uspny//2,uspnx//2:]
    imgd=np.diagonal(uspmap[:,:uspnx//2,:uspny//2], axis1=1, axis2=2)
    
    #plt.scatter(x, EigenVals_1, c=OrbVals_1, lw=0, s=10,zorder=1)

    axs[0][1].imshow(imgd,extent=[-0.707, 0.0, uspbiaslower, uspbiasupper],cmap='Greys',origin='lower',vmin=0,vmax=np.max(imgh))
    axs[0][1].imshow(imgh,extent=[0, 0.5, uspbiaslower, uspbiasupper],cmap='Greys',origin='lower',vmin=0,vmax=np.max(imgh))
    #axs.set_rasterized(True)
    
 
    axs[0][1].set_xticks([-0.707,0,0.5],['X',r'$\Gamma$','M'])
    axs[0][1].set_yticks([-0.4,-0.2,0,0.2,0.4])
    axs[0][1].set_xlim([-0.707,0.5])
    axs[0][1].set_ylim([uspbiaslower,uspbiasupper])
    axs[0][1].axhline(0,c='k',lw=0.5)
    #axs[0][1].axvline(0,c='k',lw=0.5)
    axs[0][1].set_xlabel(r'$\mathbf{k}$')
    axs[0][1].set_ylabel('Energy (eV)')
    
    qpimap,qpimapheader=Read_idl(qpifftmapname)
    qpibiaslower=float(qpimapheader[9])
    qpibiasupper=float(qpimapheader[10])
    #axs.set_rasterized(True)
    qpinx = int(qpimapheader[2]) #LDOS pixels x
    qpiny = int(qpimapheader[3]) #LDOS pixels y
    qpilayers=int(qpimapheader[4])
    fermilayer=int((0.0-qpibiaslower)/(qpibiasupper-qpibiaslower)*qpilayers)
    axs[1][0].imshow(qpimap[fermilayer,:,:],extent=[-2, 2, -2, 2],cmap='Greys',origin='lower',vmin=0,vmax=0.0003*np.max(qpimap[fermilayer,:,:]))

    axs[1][0].scatter([0],[0],s=2,c='k')
    axs[1][0].text(0.0, -0.08, '(0,0)',fontsize=12, va='top',ha='center')
    axs[1][0].text(0.0, 1.04, '(0,1)',fontsize=12, va='bottom',ha='center')
    axs[1][0].text(1.00, 1.04, '(1,1)',fontsize=12, va='bottom',ha='center')
    axs[1][0].plot([0,0],[0,1],c='r')
    axs[1][0].plot([0,1],[0,1],c='r')
    #axs[1][1].axhline(0,c='k',lw=0.5)
    #axs[1][1].axvline(0,c='k',lw=0.5)
    axs[1][0].set_xticks([-1.0,0,1.0])
    axs[1][0].set_yticks([-1.0,0,1.0])
  
    axs[1][0].set_xlim([-1.5,1.5])
    axs[1][0].set_ylim([-1.5,1.5])
    axs[1][0].set_xlabel(r'$q_x$ ($2\pi/a$)')
    axs[1][0].set_ylabel(r'$q_y$ ($2\pi/a$)')
    qpicut=qpimap[:,qpinx//2,qpiny//2:3*qpiny//4]
    qpidiag=np.diagonal(qpimap[:,qpinx//4:qpinx//2,qpiny//4:qpiny//2], axis1=1, axis2=2)
    
    im2=axs[1][1].imshow(qpicut,extent=[0, 1, qpibiaslower, qpibiasupper],aspect='auto',cmap='Greys',origin='lower',vmin=0,vmax=0.0003*np.max(qpicut))
    axs[1][1].imshow(qpidiag,extent=[-1.414, 0, qpibiaslower, qpibiasupper],aspect='auto',cmap='Greys',origin='lower',vmin=0,vmax=0.0003*np.max(qpicut))
    axs[1][1].axhline(0,c='k',lw=0.5)
    #axs[1][2].axvline(0,c='k',lw=0.5)
    axs[1][1].set_xticks([-1.414,0.0,1.0],['(1,1)','0','(0,1)'])
    axs[1][1].set_yticks([-0.4,-0.2,0,0.2,0.4])
  
    axs[1][1].set_xlim([-1.414,1.0])
    axs[1][1].set_ylim([-0.5,0.5])

    axs[1][1].set_xlabel(r'$\mathbf{q}$')
    axs[1][1].set_ylabel('Bias (V)')
 
    cbar=plt.colorbar(im2, ax=axs[1][1],ticks=[0,0.05*np.max(qpicut)])
    cbar.ax.set_yticklabels(['0',''])
    #axs[j][i].axhline(0, 0, count, lw=2, color='black')

    #filename = qpihorizname[:-4]+".pdf"
    #plt.savefig(filename, dpi=300)
    #plt.savefig(filename, bbox_inches='tight',transparent=True)
    #plt.show()


    ###############################
    #Plotting the main dispersions - diag #
    ###############################
    #fig=plt.figure(figsize=(5,5))
    #axs=fig.add_subplot(111)
    #uspmaprotated=scipy.ndimage.rotate(uspmap, angle=45, axes=(2,1))
    
    #sz,sx, sy = uspmaprotated.shape
    #imgd=uspmaprotated[:,sx//2,sy//2:3*sy//4]   
    
    #plt.scatter(x, EigenVals_1, c=OrbVals_1, lw=0, s=10,zorder=1)
 
    #axs[0][0].imshow(imgd,extent=[0, 1, biaslower, biasupper],aspect='auto',origin='lower',cmap='Greys',vmin=0,vmax=0.5*np.max(imgd))
    #axs.set_rasterized(True)
 
    #axs[0][0].set_xticks([0,1])
    #axs[0][0].set_yticks([-0.4,-0.2,0,0.2,0.4])
  
    #for i in range(1,6):
    #    plt.vlines(x=kpts*i, ymin=Ymin, ymax=Ymax)
    #plt.vlines(x=3 * self.kpoints, ymin=-2, ymax=2)
                #kpts=indarray[-1]
    #axs[0][0].set_xlim([0.0,1.0])
    #axs[0][0].set_ylim([-0.5,0.5])

    #axs[0][0].set_ylabel(r'$E (\mathrm{meV})$')
    
    #if(cutgxname is not None):
    #    axs[2][0].imshow(imgd,extent=[0, 1, biaslower, biasupper],aspect='auto',origin='lower',cmap='Greys',vmin=0,vmax=0.5*np.max(imgd))
        #axs.set_rasterized(True)
 
    #    axs[2][0].set_xticks([0,1])
    #    axs[2][0].set_yticks([-0.4,-0.2,0,0.2,0.4])
  
        #for i in range(1,6):
        #    plt.vlines(x=kpts*i, ymin=Ymin, ymax=Ymax)
        #plt.vlines(x=3 * self.kpoints, ymin=-2, ymax=2)
        #kpts=indarray[-1]
     #   axs[2][0].set_xlim([0.0,1.0])
     #   axs[2][0].set_ylim([-0.5,0.5])
     #   axs[2][0].set_ylabel(r'$E (\mathrm{meV})$')
     #   indarray,energies,store,spinstore=loadcut(cutgxname)
     #   indarray=indarray/indarray[-1]
     #   axs[2][0].scatter(indarray, energies, c=store, lw=0, s=10,zorder=1)
     #   axs[1][0].scatter(indarray, energies, c=store, lw=0, s=10,zorder=1)
    
    #plt.xticks([0, kpts,2*kpts], [r'$\overline{\Gamma}$', r'$\overline{\mathrm{M}}$', r'$\overline{\mathrm{X}}$'])
    
    #axs[0][0].set_xlabel(r'$\mathbf{k}$')
    #axs[j][i].axhline(0, 0, count, lw=2, color='black')
                
    #plt.scatter(x, EigenVals_1, c=OrbVals_1, lw=0, s=10,zorder=1)
 
    #axs.set_rasterized(True)
 
    #axs[1][0].set_xticks([0,1])
    #axs[1][0].set_yticks([-0.4,-0.2,0,0.2,0.4])
  
    #for i in range(1,6):
    #    plt.vlines(x=kpts*i, ymin=Ymin, ymax=Ymax)
    #plt.vlines(x=3 * self.kpoints, ymin=-2, ymax=2)
                #kpts=indarray[-1]
    #axs[1][0].set_xlim([0.0,1.0])
    #axs[1][0].set_ylim([-0.5,0.5])
    #axs[rows-1][0].set_xlabel(r'$\mathbf{k}$')
    #axs[1][0].set_ylabel(r'$E (\mathrm{meV})$')
    #plt.xticks([0, kpts,2*kpts], [r'$\overline{\Gamma}$', r'$\overline{\mathrm{M}}$', r'$\overline{\mathrm{X}}$'])
               

    #filename = qpidiagname[:-4]+".pdf"
    #plt.savefig(filename, dpi=300)
    #plt.savefig(filename, bbox_inches='tight',transparent=True)
    #plt.show()

    #plt.savefig(filename, dpi=300)
    for i in range(len(axs)):
        for j in range(len(axs[i])):
            label=f'({chr(i*len(axs[i])+j+97)})'
            axs[i][j].text(
                0.0, 1.0, label, transform=(
                    axs[i][j].transAxes + ScaledTranslation(-64/72, -8/72, fig.dpi_scale_trans)),
                fontsize=18, va='bottom')
    plt.savefig('1nnzfigure.pdf', bbox_inches='tight',transparent=True)

def main():
    ###To run python3 seedname Fermi_Energy Ymin Ymax
    ### e.g python3 Sr2RuO4_hr.dat 0.0 -1 1
    #cutname = str(sys.argv[1]) #Name of Hamiltonian file e.g Sr2RuO4.dat
    MkPlot('spf-bulk.idl','surfacefft.idl')
        
main()
