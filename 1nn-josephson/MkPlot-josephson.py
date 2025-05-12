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

def MkPlot(jmapname,qpimapname):
    
    font = {'family' : 'sans-serif',
        'sans-serif' : 'Arial',
        'size'   : 16}
    plt.rc('font', **font)

    ###############################
    #Plotting the main dispersions - horiz#
    ###############################
    fig=plt.figure(figsize=(7.5,5))
    gs=fig.add_gridspec(2,2,hspace=0.35)
    #gs=fig.add_gridspec(2,3, hspace=0,wspace=0)

    axs=gs.subplots()
    #(sharex='col', sharey='row')
    axs[0, 0].axis('off')
    qpimap,qpimapheader=Read_idl(qpimapname)
    qpibiaslower=float(qpimapheader[9])*1000.0
    qpibiasupper=float(qpimapheader[10])*1000.0
    jmap,jmapheader=Read_idl(jmapname)
    #axs.set_rasterized(True)
    qpinx = int(qpimapheader[2]) #LDOS pixels x
    qpiny = int(qpimapheader[3]) #LDOS pixels y
    qpilayers=int(qpimapheader[4])

    icmap=np.sum(jmap,axis=0)
 
    #fermilayer=int((-qpibiaslower)/(qpibiasupper-qpibiaslower)*qpilayers)
    vmin1=np.min(qpimap[-1,:,:])
    vmax1=np.max(qpimap[-1,:,:])
    im1=axs[0][1].imshow(qpimap[-1,:,:],extent=[-1, 1, -1, 1],cmap='Grays_r',origin='lower',vmin=vmin1,vmax=vmax1)
    cbar1=plt.colorbar(im1, ax=axs[0][1],ticks=[vmin1,vmax1])
    cbar1.ax.set_yticklabels(['min','max'])
    cbar1.ax.set_ylabel(r'$\rho$')
    axs[0][1].set_xticks([])
    axs[0][1].set_yticks([])
    axs[0][1].set_xlabel(r'$r_x$')
    axs[0][1].set_ylabel(r'$r_y$')
    axs[0][1].set_xlim([-1.0,1.0])
    axs[0][1].set_ylim([-1.0,1.0])
    for i in range(4):
        for j in range(4):
            axs[0][1].scatter([-i*0.25+1.0/128.0],[-j*0.25+1.0/128.0],c='k',s=16,linewidth=1,marker='+')
    vmin=np.min(icmap)
    vmax=np.max(icmap)
    im2=axs[1][0].imshow(icmap[:,:],extent=[-1, 1, -1, 1],cmap='seismic_r',origin='lower',vmin=vmin,vmax=vmax)
    cbar=plt.colorbar(im2, ax=axs[1][0],ticks=[vmin,vmax])
    cbar.ax.set_yticklabels(['-','+'])
    cbar.ax.set_ylabel(r'$I_\mathrm{c}$')
    axs[1][0].set_xticks([])
    axs[1][0].set_yticks([])
  
    axs[1][0].set_xlim([-1.0,1.0])
    axs[1][0].set_ylim([-1.0,1.0])
    axs[1][0].set_xlabel(r'$r_x$')
    axs[1][0].set_ylabel(r'$r_y$')
    for i in range(4):
        for j in range(4):
            axs[1][0].scatter([-i*0.25+1./128.],[-j*0.25+1./128.],c='k',s=16,linewidth=1,marker='+')
    icmap=icmap*icmap
    vmin=0
    vmax=np.max(icmap)
    im3=axs[1][1].imshow(icmap[:,:],extent=[-1, 1, -1, 1],cmap='Blues',origin='lower',vmin=vmin,vmax=vmax)
    
    cbar2=plt.colorbar(im3, ax=axs[1][1],ticks=[vmin,vmax])
    cbar2.ax.set_yticklabels([0,'max'])
    cbar2.ax.set_ylabel(r'$I_\mathrm{c}^2$')
    axs[1][1].set_xticks([])
    axs[1][1].set_yticks([])
  
    axs[1][1].set_xlim([-1.0,1.0])
    axs[1][1].set_ylim([-1.0,1.0])
    axs[1][1].set_xlabel(r'$r_x$')
    axs[1][1].set_ylabel(r'$r_y$')
    
    
    for i in range(4):
        for j in range(4):
            axs[1][1].scatter([-i*0.25+1.0/128.],[-j*0.25+1.0/128.],c='k',s=16,linewidth=1,marker='+')
    
    
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
                    axs[i][j].transAxes + ScaledTranslation(-56/72, -8/72, fig.dpi_scale_trans)),
                fontsize=18, va='bottom', fontfamily='Arial')
    plt.savefig('1nn-josephson.pdf', bbox_inches='tight',transparent=True)

def main():
    ###To run python3 seedname Fermi_Energy Ymin Ymax
    ### e.g python3 Sr2RuO4_hr.dat 0.0 -1 1
    #cutname = str(sys.argv[1]) #Name of Hamiltonian file e.g Sr2RuO4.dat
    MkPlot('1nn-dx2j.idl','1nn-dx2jn.idl')
        
main()
