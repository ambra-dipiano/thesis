# MIT License
# Copyright (c) 2019, 2020 Ambra Di Piano
# ---------------------------------------
# =====================
# !!! SPECTRA DN LC !!!
# =====================

import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import seaborn as sns

pathin = '/home/ambra/Desktop/cluster-morgana/run0406_test/run0406/' 
path = pathin + 'run0406_ID000126/'
png = pathin + 'png/'
template = 'run0406_ID000126_ebl.fits'

# load data ---!
hdul =  fits.open(pathin + template)
print(hdul.info())

# energybins [GeV]
energy=np.array(hdul[1].data)
# timebins [s]
time=np.array(hdul[2].data)
# spectra [fotoni/GeV/cm^2/s]
spectra=np.array(hdul[3].data)
# ebl [fotoni/GeV/cm^2/s]
ebl=np.array(hdul[4].data)

hdul.close()

Nt=len(time)
Ne=len(energy)

# TIME GRID ---!
t=[0.0 for x in range(Nt+1)]
for i in range(Nt-1):
    t[i+1]=time[i][0]+(time[i+1][0]-time[i][0])/2
# last bin ---!
t[Nt]=time[Nt-1][0]+(time[Nt-1][0]-t[Nt-1])

# ENERGY GRID ---!
en=[1.0 for x in range(Ne+1)]
for i in range(Ne-1):
    en[i+1]=energy[i][0]+(energy[i+1][0]-energy[i][0])/2
# last bin ---!
en[Ne]=energy[Ne-1][0]+(energy[Ne-1][0]-en[Ne-1])

# energy ranges
LST = (min(en, key=lambda x:abs(x-30)), min(en, key=lambda x:abs(x-150)))
MST = (min(en, key=lambda x:abs(x-150)), min(en, key=lambda x:abs(x-10000)))
SST = (min(en, key=lambda x:abs(x-1000)), min(en, key=lambda x:abs(x-300000)))
CTA = (min(en, key=lambda x:abs(x-30)), min(en, key=lambda x:abs(x-150000)))
below = (min(en, key=lambda x:abs(x-1)), min(en, key=lambda x:abs(x-30)))

# time bins ---!
medium = min(en, key=lambda x:abs(x-30))
worst = min(en, key=lambda x:abs(x-50))
hundred = min(en, key=lambda x:abs(x-100))
thousand = min(en, key=lambda x:abs(x-1000))
tbins = [medium, worst, hundred, thousand]

# =====================
# !!! PLOT TEMPLATE !!!
# =====================

sns.set()

# FLUX SPECTRA ---!
f=[]
f8=[]
f9=[]
f10=[]
for i in range(Nt):
    f.append(0.0)
    f8.append(0.0)
    f9.append(0.0)
    f10.append(0.0)
    for j in range(Ne):
        f[i]=f[i]+spectra[i][j]*(en[j+1]-en[j])
        if en[j] <= SST[1] and en[j] >= SST[0]:
            f8[i]=f8[i]+spectra[i][j]*(en[j+1]-en[j])
        if en[j] <= CTA[1] and en[j] >= CTA[0]:
            f9[i]=f9[i]+spectra[i][j]*(en[j+1]-en[j])
        if en[j] <= MST[1] and en[j] >= MST[0]:
            f10[i]=f10[i]+spectra[i][j]*(en[j+1]-en[j])

# FLUX EBL ---!
f2=[]
f3=[]
f4=[]
f5=[]
f6=[]
f7=[]
for i in range(Nt):
    f2.append(0.0)
    f3.append(0.0)
    f4.append(0.0)
    f5.append(0.0)
    f6.append(0.0)
    f7.append(0.0)
    for j in range(Ne):
        f2[i]=f2[i]+ebl[i][j]*(en[j+1]-en[j])
        if en[j] <= LST[1] and en[j] >= LST[0]:
            f3[i]=f3[i]+ebl[i][j]*(en[j+1]-en[j])
        if en[j] <= MST[1] and en[j] >= MST[0]:
            f4[i]=f4[i]+ebl[i][j]*(en[j+1]-en[j])
        if en[j] <= CTA[1] and en[j] >= CTA[0]:
            f5[i]=f5[i]+ebl[i][j]*(en[j+1]-en[j])
        if en[j] < LST[0]:
            f6[i]=f6[i]+ebl[i][j]*(en[j+1]-en[j])
        if en[j] <= SST[1] and en[j] >= SST[0]:
            f7[i]=f7[i]+ebl[i][j]*(en[j+1]-en[j])
          
# some t spectra ---!
fig1=plt.figure(3)
t1=27
for i in range(2):
    x=[]
    y=[]
    z=[]
    for j in range(Ne):
        x.append(energy[j][0])
        y.append(spectra[t1+i*20][j])
        z.append(ebl[t1+i*20][j])
    plt.loglog(x, y, '-', label='%0.2fs-%0.2fs' %(time[t1+i*20][0], time[t1+i*20+1][0]))
    plt.loglog(x, z, '-.', label='" " EBL adsorbed')
    plt.title('template spectra')
    plt.xlabel('E (GeV)')
    plt.ylabel('d$\Phi$/dE (ph/cm$^2$/s/GeV)')
plt.axvline(30, c='k', ls='--')
plt.text(40, 1e-8, '30 GeV', color='k')
plt.legend(loc=3) #if i == 0 else None
#plt.xlim([0.1, 1e4])
plt.tight_layout()
fig1.savefig(png + 'template_spectra_ebl.png')
plt.show()
  
# lightcurves ---!
fig2=plt.figure(2)
ax = plt.subplot(111, yscale='log')
#plt.plot(time, f, label='no EBL')
plt.plot(time, f9, label='CTA (30GeV-150TeV)')
plt.plot(time, f5, ls='-.', label='CTA+EBL')
plt.plot(time, f10, label='MST (150GeV-10TeV)')
plt.plot(time, f4, ls='-.', label='MST+EBL')
plt.title('template lightcurve')
plt.xlabel('time (s)')
plt.ylabel('F (ph/cm2/s)')
plt.legend()
plt.ylim([1e-11, 1e-8])
plt.xlim([-0.5, 2e3])
plt.tight_layout()
fig2.savefig(png + 'template_lightcurve.png')
plt.show()


