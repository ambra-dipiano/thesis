import gammalib
import ctools
import cscripts
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import numpy as np


hdul =  fits.open('/home/ambra/Desktop/cluster-morgana/run0406_lc/run0406/run0406_ID000126.fits')
#hdul.info()

## 41 energybins [GeV]
energy=np.array(hdul[1].data)
# 71 timebins [s]
time=np.array(hdul[2].data)
# 71 spettri  [array 41 flussi x 71 timebins, ph/cm2sGeV]
spectra=np.array(hdul[3].data)

Nt=len(time)
Ne=len(energy)
print('!!! check ', time)

#plot flux at 100 GeV vs time
x500=[]
y500=[]
x100=[]
y100=[]
x30=[]
y30=[]
c = 0
for i in range(Nt):
    if time[i][0] <= 150:
        c += 1
        print(i, time[i][0],spectra[i][0])
        x500.append(time[i][0])
        y500.append(spectra[i][27])
        x100.append(time[i][0])
        y100.append(spectra[i][20])
        x30.append(time[i][0])
        y30.append(spectra[i][15])

bin_stop = c
plt.loglog(x500,y500,'b',label='500 GeV')
plt.loglog(x100,y100,'g',label='100 GeV')
plt.loglog(x30,y30,'r',label='30 GeV')
plt.xlabel('Time[s]')
plt.ylabel('Photon Flux [ph/cm2sGeV]')
plt.legend()
plt.show()

# DEFINISCO LA GRIGLIA TEMPORALE
t=[0.0 for x in range(Nt+1)]
for i in range(Nt-1):
    t[i+1]=time[i][0]+(time[i+1][0]-time[i][0])/2
# definisco il tmax dell'ultimo bin
t[Nt]=time[Nt-1][0]+(time[Nt-1][0]-t[Nt-1])


# DEFINISCO LA GRIGLIA DI ENERGIE
en=[1.0 for x in range(Ne+1)]
for i in range(Ne-1):
    en[i+1]=energy[i][0]+(energy[i+1][0]-energy[i][0])/2
# definisco il tmax dell'ultimo bin
en[Ne]=energy[Ne-1][0]+(energy[Ne-1][0]-en[Ne-1])


# CALCOLO IL FLUSSO IN CIASCUN BIN temporale
f=[0.0 for x in range(Nt)]
for i in range(Nt):
    f[i]=0
    for j in range(Ne):
        #print energy[j][0],spectra[i][j]
        f[i]=f[i]+spectra[i][j]*(en[j+1]-en[j])


# Plot degli spettri nei primi 10 bin temporali:
for i in range(10):
    x=[]
    y=[]
    for j in range(Ne):
        x.append(energy[j][0])
        y.append(spectra[i][j])
    f1=plt.figure(1)
    plt.loglog(x,y)
    plt.title('run0406_ID000126')
    plt.xlabel('E[GeV]')
    plt.ylabel('F(E)[ph/cm2sGeV] during first 10 time bins')
    plt.savefig('template_spectra.png')
    f1.show()


# PLOTTA LA CURVA DI LUCE
f2=plt.figure(2)
plt.loglog(time,f)
plt.title('run0406_ID000126')
plt.xlabel('time[s]')
plt.ylabel('F[ph/cm2s]')
plt.savefig('template_lightcurve.png')
f2.show()

#raw_input()
