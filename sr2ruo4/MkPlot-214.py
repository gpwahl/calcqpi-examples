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

def swaplayers(data,header):
    newdata=np.flip(np.copy(data),axis=0)
    newheader=header.copy()
    newheader[9],newheader[10]=header[10],header[9]
    return newdata,newheader


def calcprfft(qpimapheader,qpimap):
    qpibiaslower=float(qpimapheader[9])
    qpibiasupper=float(qpimapheader[10])
    qpilayers=int(qpimapheader[4])
    zerobias=int((-qpibiaslower)/(qpibiasupper-qpibiaslower)*qpilayers)
    for i in range(len(qpimap)):
       qpimap[i]=qpimap[i]-np.average(qpimap[i])
    cqpimap=np.fft.fftshift(np.fft.fft2(qpimap),axes=(1,2))
    qpimap=np.abs(cqpimap)
    for ref in range(zerobias):
        i=len(qpimap)-ref-1
        qpimap[i]=np.real(cqpimap[i,:,:]*np.abs(cqpimap[ref,:,:])/cqpimap[ref,:,:])
    return qpimap

def calcfft(qpimap):
    for i in range(len(qpimap)):
       qpimap[i]=qpimap[i]-np.average(qpimap[i])
    qpimap=np.fft.fftshift(np.abs(np.fft.fft2(qpimap)),axes=(1,2))
    return qpimap

def MkPlot(qpi01mapname,qpi10mapname):
    
    font = {'family' : 'sans-serif',
        'sans-serif' : 'Arial',
        'size'   : 16}
    plt.rc('font', **font)

    ###############################
    #Plotting the main dispersions - horiz#
    ###############################
    fig=plt.figure(figsize=(14.5,4))
    gs=fig.add_gridspec(1,3,wspace=0.3)
    #gs=fig.add_gridspec(2,3, hspace=0,wspace=0)

    axs=gs.subplots()
    #(sharex='col', sharey='row')
     
    qpi01map,qpi01mapheader=Read_idl(qpi01mapname)
    qpi10map,qpi10mapheader=Read_idl(qpi10mapname)
    qpibiaslower=float(qpi01mapheader[9])*1000.0
    qpibiasupper=float(qpi01mapheader[10])*1000.0
    #axs.set_rasterized(True)
    qpinx = int(qpi01mapheader[2]) #LDOS pixels x
    qpiny = int(qpi01mapheader[3]) #LDOS pixels y
    qpilayers=int(qpi01mapheader[4])

    avespec=np.average(qpi01map,axis=(1,2))+np.average(qpi10map,axis=(1,2))
    avespec=avespec/avespec[-1]
    bias=np.linspace(qpibiaslower,qpibiasupper,qpilayers,endpoint=True)
    axs[0].scatter(bias,avespec,c='k')
    axs[0].set_xticks([-8,-4,0,4,8])
    axs[0].set_xlim([-10,10])
    axs[0].set_xlabel('Bias (mV)')
    axs[0].set_ylabel(r'$\rho(eV)$ (rel. u.)')
    refenergy=-2.8
    axs[0].axvline(refenergy,c='b',linestyle='dashed')
  
    fermilayer=int((refenergy-qpibiaslower)/(qpibiasupper-qpibiaslower)*qpilayers)
    fftqpimap=calcfft(qpi01map)+calcfft(qpi10map)
    vmax=0.01*np.max(fftqpimap[fermilayer,:,:])
    vmin=-vmax
    axs[1].imshow(fftqpimap[fermilayer,:,:],extent=[-2, 2, -2, 2],cmap='Blues',origin='lower',vmin=0.0,vmax=vmax)
    axs[1].set_xticks([-1.0,0,1.0])
    axs[1].set_yticks([-1.0,0,1.0])
    axs[1].set_xlabel(r'$q_x (2\pi/a)$')
    axs[1].set_ylabel(r'$q_y (2\pi/a)$')
    axs[1].set_xlim([-1.5,1.5])
    axs[1].set_ylim([-1.5,1.5])
    axs[1].plot([0,0],[0,1],color='grey',linestyle='dashed')
    axs[1].plot([0,1],[0,1],color='grey',linestyle='dashed')
    qpicut=fftqpimap[:,qpinx//2,qpiny//2:3*qpiny//4]
    qpidiag=np.diagonal(fftqpimap[:,qpinx//4:qpinx//2,qpiny//4:qpiny//2], axis1=1, axis2=2)
    
    axs[2].imshow(qpicut,extent=[0, 1, qpibiaslower, qpibiasupper],aspect='auto',cmap='Blues',origin='lower',vmin=0.0,vmax=vmax)
    im2=axs[2].imshow(qpidiag,extent=[-1.414, 0, qpibiaslower, qpibiasupper],aspect='auto',cmap='Blues',origin='lower',vmin=0.0,vmax=vmax)
    axs[2].set_xticks([-1.414,0.0,1.0],['(1,1)','0','(1,0)'])
    axs[2].set_yticks([-20,0,20])
  
    axs[2].set_xlim([-1.414,1.0])
    axs[2].set_ylim([qpibiaslower,qpibiasupper])

    axs[2].set_xlabel(r'$\mathbf{q} (2\pi/a)$')
    axs[2].set_ylabel('Bias (mV)')
    
    cbar=plt.colorbar(im2, ax=axs[2],ticks=[0,vmax])
    cbar.ax.set_yticklabels(['0','+'])
    
    for j in range(len(axs)):
        label=f'({chr(j+97)})'
        axs[j].text(
            0.0, 1.0, label, transform=(
                axs[j].transAxes + ScaledTranslation(-56/72, -8/72, fig.dpi_scale_trans)),
            fontsize=18, va='bottom', fontfamily='Arial')
    plt.savefig('214figure.pdf', bbox_inches='tight',transparent=True)

def main():
    ###To run python3 seedname Fermi_Energy Ymin Ymax
    ### e.g python3 Sr2RuO4_hr.dat 0.0 -1 1
    #cutname = str(sys.argv[1]) #Name of Hamiltonian file e.g Sr2RuO4.dat
    MkPlot('214qpi01.idl','214qpi10.idl')
        
main()
