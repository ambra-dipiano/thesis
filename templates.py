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
erange = 'fullsys'

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

# energy in erg (1 GeV = 0.00160218 erg) ---!
erg = np.array(en)/0.00160218

# energy ranges
if erange == 'fullsys':
    LST = (min(en, key=lambda x:abs(x-20)), min(en, key=lambda x:abs(x-150)))
    lLST = '20GeV-150GeV'
    MST = (min(en, key=lambda x:abs(x-150)), min(en, key=lambda x:abs(x-5000)))
    lMST = '150GeV-5TeV'
    SST = (min(en, key=lambda x:abs(x-5000)), min(en, key=lambda x:abs(x-300000)))
    lSST = '5TeV-300TeV'
else:
    LST = (min(en, key=lambda x:abs(x-20)), min(en, key=lambda x:abs(x-3000)))
    lLST = '20GeV-3TeV'
    MST = (min(en, key=lambda x:abs(x-80)), min(en, key=lambda x:abs(x-50000)))
    lMST = '80GeV-50TeV'
    SST = (min(en, key=lambda x:abs(x-1000)), min(en, key=lambda x:abs(x-300000)))
    lSST = '1TeV-300TeV'
CTA = (min(en, key=lambda x:abs(x-30)), min(en, key=lambda x:abs(x-150000)))
lCTA = '30GeV-150TeV'
MAGIC = (min(en, key=lambda x:abs(x-300)), min(en, key=lambda x:abs(x-1000)))
lMAGIC = '30GeV-1TeV'
MAGIC2 = (min(en, key=lambda x: abs(x - 300)), min(en, key=lambda x: abs(x - 10000)))
lMAGIC2 = '30GeV-1TeV'
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
f11=[]
for i in range(Nt):
    f.append(0.0)
    f8.append(0.0)
    f9.append(0.0)
    f10.append(0.0)
    f11.append(0.0)
    for j in range(Ne):
        f[i]=f[i]+spectra[i][j]*(en[j+1]-en[j])
        if en[j] <= SST[1] and en[j] >= SST[0]:
            f8[i]=f8[i]+spectra[i][j]*(en[j+1]-en[j])
        if en[j] <= CTA[1] and en[j] >= CTA[0]:
            f9[i]=f9[i]+spectra[i][j]*(en[j+1]-en[j])
        if en[j] <= MST[1] and en[j] >= MST[0]:
            f10[i] = f10[i] + spectra[i][j] * (en[j + 1] - en[j])
        if en[j] <= LST[1] and en[j] >= LST[0]:
            f11[i]=f11[i]+spectra[i][j]*(en[j+1]-en[j])

# FLUX EBL ---!
f2=[]
f3=[]
f4=[]
f5=[]
f6=[]
f7=[]
ferg=[]
ferg2=[]
f12=[]
f13=[]
for i in range(Nt):
    f2.append(0.0)
    f3.append(0.0)
    f4.append(0.0)
    f5.append(0.0)
    f6.append(0.0)
    f7.append(0.0)
    ferg.append(0.0)
    ferg2.append(0.0)
    f12.append(0.0)
    f13.append(0.0)
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
        if en[j] <= CTA[1] and en[j] >= CTA[0]:
            ferg[i]=ferg[i]+ebl[i][j]*0.0016022*(erg[j+1]-erg[j])
        if en[j] <= MAGIC[1] and en[j] >= MAGIC[0]:
            ferg2[i]=ferg2[i]+ebl[i][j]*0.0016022*(erg[j+1]-erg[j])
        if en[j] <= MAGIC[1] and en[j] >= MAGIC[0]:
            f12[i]=f12[i]+ebl[i][j]*(en[j+1]-en[j])
        if en[j] <= MAGIC2[1] and en[j] >= MAGIC2[0]:
            f13[i]=f13[i]+ebl[i][j]*(en[j+1]-en[j])

# some t spectra ---!
fig1=plt.figure(3)
t1=27
for i in range(3):
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
fig1.savefig(png + 'template_spectra_ebl_3.png')
plt.show()
  
# lightcurves ---!
fig2=plt.figure(2)
ax = plt.subplot(111, yscale='log')
#plt.plot(time, f, label='no EBL')
plt.plot(time, f9, label='CTA (%s)' %lCTA)
plt.plot(time, f5, ls='-.', label='CTA+EBL')
plt.plot(time, f11, label='LST (%s)' %lLST)
plt.plot(time, f3, ls='-.', label='LST+EBL')
plt.plot(time, f10, label='MST (%s)' %lMST)
plt.plot(time, f4, ls='-.', label='MST+EBL')
plt.plot(time, f8, label='SST (%s)' %lSST)
plt.plot(time, f7, ls='-.', label='SST+EBL')
if erange == 'core':
    ranges = '(full system sensitivity range)'
    suffix = 'fullsys'
else:
    ranges = '(required range)'
    suffix = 'req'
plt.title('template lightcurve '+ranges)
plt.xlabel('time (s)')
plt.ylabel('F (ph/cm2/s)')
plt.legend()
plt.ylim([1e-13, 1e-8])
plt.xlim([-0.5, 15e2])
plt.tight_layout()
fig2.savefig(png + 'template_lightcurve_%s.png' %suffix)
plt.show()

# lightcurves [erg] ---!
fig2=plt.figure(2)
ax = plt.subplot(111, yscale='log', xscale='linear')
#plt.plot(time, f, label='no EBL')
#plt.plot(time, f2, label='CTA (%s)' %lCTA)
plt.plot(time, f12, label='EBL absorbed (%s)' %lMAGIC)
#plt.plot(time, f13, label='EBL absorbed (%s)' %lMAGIC2)
ranges = '(E range of GRB190114C obs. by MAGIC)'
suffix = 'MAGIC'
plt.title('template lightcurves'+ranges)
plt.xlabel('time (s)')
plt.ylabel('F (ph/cm2/s)')
plt.legend()
plt.ylim([1e-11, 1e-9])
plt.xlim([0, 200])
plt.tight_layout()
fig2.savefig(png + 'template_lightcurve_%s.png' %suffix)
plt.show()
