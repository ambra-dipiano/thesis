# ============================ #
# MODULE OF PLOTTING FUNCTIONS #
# ============================ #

# IMPORTS ---!
import matplotlib.pyplot as plt
import pyregion
import seaborn as sns
from astropy.io import fits
from matplotlib.colors import SymLogNorm
from matplotlib import rc
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Rectangle
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename


# WIP ---!
def handleReg(reg, col='black') :
  '''
  Changes attributes of region
  :param reg: filename (str)
  :param col: color (str)
  :return:
  '''
  r = pyregion.open(reg)
  r[0].attr[1]['color'] = col
  #print(r[0].attr[1]['color'])

  # NEED TO SAVE CHANGES ---!

  return r

def showSkymap(skymap, reg='none', col='black', suffix='none', title='none', show=True):
  ''''
  :param:
  skymap = skymap fits file (str)
  index = index of the extention storing the data (int)
  reg = region reg file (str)
  col = color to plot the region with (str)
  suffix = suffix added to image filename (str)
  show = if calling image show (bool)
  :return:
  '''

  with fits.open(skymap) as hdul:
    wcs = WCS(hdul[0].header)
    data = hdul[0].data

  # PLOT ---!
  #plt.rc('text', usetex=True)
  ax = plt.subplot(111)

  # handle region ---!
  if reg != 'none' :
    r = pyregion.open(reg).as_imagecoord(hdul[0].header)
    for i in range(len(r)) :
      r[i].attr[1]['color'] = col
      patch_list, text_list = r.get_mpl_patches_texts()
      for p in patch_list:
        ax.add_patch(p)
      for t in text_list:
        ax.add_artist(t)


  plt.subplot(projection=wcs)
  plt.imshow(data, cmap='jet', norm=SymLogNorm(1), interpolation='gaussian')
  plt.grid(color='white', ls='solid')
  plt.xlabel('R.A. (deg)')
  plt.ylabel('Dec (deg)')
  plt.title('skymap')  if title == 'none' else plt.title(title)
  plt.colorbar().set_label('cts')

  # save fig ---!
  if suffix != 'none' :
    plt.savefig(skymap.replace('.fits', '_%s.png' % suffix))
  else :
    plt.savefig(skymap.replace('.fits', '.png'))

  # show fig ---!
  plt.show() if show == True else None
  plt.close()

  return

def showResmap(resmap, reg='none', col='black', suffix='none', title='none', show=True):
  '''
  :param:
  resmap = resmap fits file (str)
  show = if calling image show (bool)
  :return:
  '''

  # info data ---!
  hdulist = fits.open(resmap)
  #print(hdulist.info())

  # get data ---!
  data = hdulist[0].data
  #print('\n shape:', data.shape, '\n type:', data.dtype.name)

  hdulist.close()

  # PLOT ---!
  ax = plt.subplot(111)
  plt.rc('text', usetex=True)

  # handle region ---!
  if reg != 'none':
    r = pyregion.open(reg).as_imagecoord(hdulist[0].header)
    for i in range(len(r)):
      r[i].attr[1]['color'] = col
      patch_list, text_list = r.get_mpl_patches_texts()
      for p in patch_list:
        ax.add_patch(p)
      for t in text_list:
        ax.add_artist(t)

  plt.imshow(data, cmap='bwr', interpolation='gaussian', filterrad=0.04)
  #plt.xlabel('R.A. (deg)')
  #plt.ylabel('Dec (deg)')
  plt.xlabel('pix')
  plt.ylabel('pix')
  plt.title('spatial residuals')  if title == 'none' else plt.title(title)
  plt.colorbar().set_label('Significance $\sigma$')

  # save fig ---!
  if suffix != 'none':
    plt.savefig(resmap.replace('.fits', '_%s.png' % suffix))
  else:
    plt.savefig(resmap.replace('.fits', '.png'))

  #show fig ---!
  plt.show() if show == True else None
  plt.close()

  return
 
def showResiduals(residuals, scaleY='log', title='none', show=True):
  '''
  :param
  residuals = residuals fits file (str)
  show = if calling image show (bool)
  :return:
  '''

  # manage hdu and info ---!
  hdulist = fits.open(residuals)
  #print(hdulist.info())

  # get data ---!
  data = hdulist[1].data
  #print('\n shape:', data.shape, '\n type:', data.dtype.name)

  # store data ---!
  Emin = data.field(0)
  Emax = data.field(1)
  cts = data.field(2)
  model = data.field(3)
  res = data.field(4)

  hdulist.close()

  # binning ---!
  en_bins = Emax - 0.5 * (Emax - Emin)

  # PLOT ---!
  plt.figure(figsize=(10,8))
  plt.rc('text', usetex=True)

  if scaleY.lower() == 'lin' :
    ax1 = plt.subplot(211, xscale='log')
  else :
    ax1 = plt.subplot(211, yscale='log', xscale='log')

  plt.plot(en_bins, cts, 'ro')
  plt.step(Emax, cts, 'r-', where='pre', label='cts')
  plt.step(Emin, model, 'g-', where='post', label='model')
  plt.xlabel('Energy (TeV)')
  plt.ylabel('cts')
  plt.title('spectral residuals')  if title == 'none' else plt.title(title)
  plt.legend(loc=0)
  plt.grid(True)

  plt.subplot(212, sharex=ax1)
  plt.plot(en_bins, res, 'b+', label='cts residuals')
  plt.axhline(y=0, c='r', lw='1', ls='--')
  plt.xlabel('Energy (TeV)')
  plt.ylabel('Residuals $\sigma$')
  plt.grid(True)

  plt.subplots_adjust(hspace=0.2)

  # save fig ---!
  plt.savefig(residuals.replace('.fits', '.png'))

  # show fig ---!
  plt.show() if show == True else None
  plt.close()

  return
 
def showButterfly(diagram, flux_pnts=0.0, fluxEn_pnts=0.0, suffix='none', title='none', show=True):
  '''
  :param
  diagram = butterfly diagram txt file (str)
  flux_pnts = list of fluxes in ph/cm^2/s/MeV (list of floats)
  fluxEn_pnts = list of energies in MeV
  show = if calling image show (bool)
  :return:
  '''
  data = np.loadtxt(diagram, delimiter=' ')

  # info ---!
  #print('\n file shape:', data.shape)

  # store data ---!
  energy = data[:, 0]
  intensity = data[:, 1]
  lerr = data[:, 2]
  uerr = data[:, 3]

  # from intensity (ph/cm^2/s/MeV) find flux in (erg/cm^2/s) ---!
  (flux, f_lerr, f_uerr) = ((intensity, lerr, uerr) * (energy ** 2)) / (6.4215 * 1e5)
  f_pnts = (np.array(flux_pnts) * (np.array(fluxEn_pnts) ** 2)) / (6.4215 * 1e5) if flux_pnts != 0.0 \
                                                                and fluxEn_pnts != 0.0 else print('plotting without flux points')

  # PLOT ---!

  ax = plt.subplot(111, yscale='log', xscale='log')
  plt.rc('text', usetex=True)

  plt.plot(energy / 1e6, flux, 'b-', alpha=1, label='best fit')
  plt.fill_between(energy / 1e6, f_lerr, f_uerr, facecolor='blue', alpha=0.3, label='errors')
  plt.scatter(np.array(fluxEn_pnts) / 1e6, f_pnts, marker='o', color='red', label='data flux pnts') if flux_pnts != 0.0 \
                                                                             and fluxEn_pnts != 0.0 else None
  plt.xlabel('Energy (TeV)')
  plt.ylabel('E $\cdot \\frac{dN}{dE}$ (erg/$cm^2$/s)')
  plt.title('butterfly diagram') if title == 'none' else plt.title(title)
  plt.legend(loc=0)
  plt.grid(True)

  # save fig ---!
  if suffix != 'none':
    plt.savefig(diagram.replace('.txt', '_%s.png' % suffix))
  else:
    plt.savefig(diagram.replace('.txt', '.png'))

  # show fig ---!
  plt.show() if show == True else None
  plt.close()

  return
 
def showSpectrum(spectrum, title='none', show=True) :
  '''
  :param spectrum: spectrum fits file (str)
  show = if calling image show (bool)
  :return:
  '''
  hdulist = fits.open(spectrum)
  #print(hdulist.info())

  # get data ---!
  data = hdulist[1].data
  #print('\n shape:', data.shape, '\n type:', data.dtype.name)

  # store data ---!
  en = data.field(0)
  en_down = data.field(1)
  en_up = data.field(2)
  flux = data.field(3)
  err_flux = data.field(4)

  hdulist.close()

  # adjusting errors ---!
  en_err = (en_up - en_down) / 2

  # PLOT ---!
  fig = plt.figure(figsize=(8,15))
  plt.rc('text', usetex=True)

  ax1 = plt.subplot(311, xscale='log')
  plt.errorbar(en, flux, yerr=err_flux, xerr=en_err, fmt='ro', label='data')
  plt.step(en, flux, 'r-', where='mid')
  plt.xlabel('Energy (TeV)')
  plt.ylabel('Flux (erg/$cm^2$/s)')
  plt.title('spectrum with errors') if title == 'none' else plt.title(title)
  plt.legend(loc=0)
  plt.grid(True)

  ax2 = plt.subplot(312, sharex=ax1)
  plt.plot(en, flux, 'ro', label='data')
  plt.step(en, flux, 'r-', where='mid')
  plt.xlabel('Energy (TeV)')
  plt.ylabel('Flux (erg/$cm^2$/s)')
  plt.title('spectrum without errors') if title == 'none' else plt.title(title)
  plt.legend(loc=0)
  plt.grid(True)

  ax3 = plt.subplot(313, sharex=ax1, yscale='log')
  plt.plot(en, flux, 'ro', label='data')
  plt.step(en, flux, 'r-', where='mid')
  plt.xlabel('Energy (TeV)')
  plt.ylabel('Flux (erg/$cm^2$/s)')
  plt.title('log spectrum without errors') if title == 'none' else plt.title(title)
  plt.legend(loc=0)
  plt.grid(True)

  plt.subplots_adjust(hspace=0.5)

  # save fig ---!
  extent = ax1.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
  fig.savefig(spectrum.replace('.fits', '_errors.png'), bbox_inches=extent.expanded(1.3, 1.3))

  extent = ax2.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
  fig.savefig(spectrum.replace('.fits', '.png'), bbox_inches=extent.expanded(1.3, 1.3))

  extent = ax3.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
  fig.savefig(spectrum.replace('.fits', '_log.png'), bbox_inches=extent.expanded(1.3, 1.3))

  # show fig ---!
  plt.show() if show == True else None
  plt.close()

  return
 
# WITH UPPER LIMITS v01 ---!
def showLightcurve(lightcurve, axisLim ='auto', title='none', show = True):
  '''
  :param lightcurve: lightcurve fits file (str)
  :param axisLim: axis range [xmin, xmax, ymin, ymax] (list)
  :return:
  '''

  hdulist = fits.open(lightcurve)
  #print(hdulist.info())

  # get data ---!
  data = hdulist[1].data
  #print('\n shape:', data.shape, '\n type:', data.dtype.name)

  # store data ---!
  t_mjd = data.field(0)  # days
  et_mjd = data.field(1)  # days
  #RA = data.field(2) # deg
  #eRA = data.field(3) # deg
  #DEC = data.field(4) # deg
  #eDEC = data.field(5) # deg
  prefact = data.field(6)  # ph/cm^2/s/MeV
  e_prefact = data.field(7)  # ph/cm^2/s/MeV
  #index = data.field(8)
  #e_index = data.field(9)
  TS = data.field(10)
  diff_uplim = data.field(11) # ph/cm^2/s/MeV
  #flux_uplim = data.field(12) # ph/cm^2/s
  #Eflux_uplim = data.field(13) # erg/cm^2/s

  hdulist.close()

  pnts = []
  e_pnts = []
  t_pnts = []
  et_pnts = []
  ul_pnts = []
  eul_pnts = []
  tul_pnts = []
  etul_pnts = []

  # list flux point or upper limit ---!
  for el in range(len(data)) :
    if TS[el] > 9 and 2.0*e_prefact[el] < prefact[el] :
    #if TS[el] > 9 :
      pnts.append(prefact[el])
      e_pnts.append(e_prefact[el])
      t_pnts.append(t_mjd[el])
      et_pnts.append(et_mjd[el])

    else :
      ul_pnts.append(diff_uplim[el])
      eul_pnts.append(0.5*diff_uplim[el])
      tul_pnts.append(t_mjd[el])
      etul_pnts.append(et_mjd[el])


  # PLOT ---!
  fig = plt.figure(figsize=(15,15))
  plt.rc('text', usetex=True)

  # figsize = (x, y)
  ax1 = plt.subplot(211, yscale='linear')
  plt.errorbar(t_pnts, pnts, xerr=et_pnts, yerr=e_pnts, fmt='ro', mec='k', label='data')
  plt.errorbar(tul_pnts, ul_pnts, xerr=[etul_pnts, etul_pnts], yerr=eul_pnts, uplims=True, fmt='ro', mec='k')
  plt.axis(axisLim) if axisLim != 'auto' else None
  plt.grid()
  plt.xlabel('t (MJD)')
  plt.ylabel('dN/dE (ph/$cm^2$/s/MeV)')
  plt.title('lightcurve') if title == 'none' else plt.title(title)
  plt.legend(loc=0)

  ax2 = plt.subplot(212, yscale='log')
  plt.errorbar(t_pnts, pnts, xerr=et_pnts, yerr=e_pnts, fmt='ro', mec='k', label='data')
  plt.errorbar(tul_pnts, ul_pnts, xerr=[etul_pnts, etul_pnts], yerr=eul_pnts, uplims=True, fmt='ro', mec='k')
  plt.axis(axisLim) if axisLim != 'auto' else None
  plt.grid()
  plt.xlabel('t (MJD)')
  plt.ylabel('dN/dE (ph/$cm^2$/s/MeV)')
  plt.title('lightcurve') if title == 'none' else plt.title(title)
  plt.legend(loc=0)

  plt.subplots_adjust(hspace=0.5)

  # save fig ---!
  extent = ax1.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
  fig.savefig(lightcurve.replace('.fits', '.png'), bbox_inches=extent.expanded(1.3, 1.3))

  extent = ax2.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
  fig.savefig(lightcurve.replace('.fits', '_log.png'), bbox_inches=extent.expanded(1.3, 1.3))

  # show fig ---!
  plt.show() if show == True else None
  plt.close()

  return
 
# WITH UPPER LIMITS v02 ---!
def showLightcurve_v02(lightcurve, axisLim ='auto', title='none', ax_scale='lin', show = True):
  '''
  :param lightcurve: lightcurve fits file (str)
  :param axisLim: axis range [xmin, xmax, ymin, ymax] (list)
  :return:
  '''

  # PLOT ---!
  fig = plt.figure(figsize=(15,15))
  plt.rc('text', usetex=True)

  # figsize = (x, y)
  ax1 = plt.subplot(211, yscale='linear') if ax_scale == 'lin' else None
  ax2 = plt.subplot(212, yscale='log') if ax_scale == 'log' else None


  for i in range(len(lightcurve)) :
    hdulist = fits.open(lightcurve)
    #print(hdulist.info())

    # get data ---!
    data = hdulist[1].data
    #print('\n shape:', data.shape, '\n type:', data.dtype.name)

    # store data ---!
    t_mjd = data.field(0)  # days
    et_mjd = data.field(1)  # days
    #RA = data.field(2) # deg
    #eRA = data.field(3) # deg
    #DEC = data.field(4) # deg
    #eDEC = data.field(5) # deg
    prefact = data.field(6)  # ph/cm^2/s/MeV
    e_prefact = data.field(7)  # ph/cm^2/s/MeV
    #index = data.field(8)
    #e_index = data.field(9)
    TS = data.field(10)
    diff_uplim = data.field(11) # ph/cm^2/s/MeV
    #flux_uplim = data.field(12) # ph/cm^2/s
    #Eflux_uplim = data.field(13) # erg/cm^2/s

    hdulist.close()

    pnts = []
    e_pnts = []
    t_pnts = []
    et_pnts = []
    ul_pnts = []
    eul_pnts = []
    tul_pnts = []
    etul_pnts = []

    # list flux point or upper limit ---!
    for el in range(len(data)) :
      if TS[el] > 9 and 2.0*e_prefact[el] < prefact[el] :
        pnts.append(prefact[el])
        e_pnts.append(e_prefact[el])
        t_pnts.append(t_mjd[el])
        et_pnts.append(et_mjd[el])

      else :
        ul_pnts.append(diff_uplim[el])
        eul_pnts.append(0.5*diff_uplim[el])
        tul_pnts.append(t_mjd[el])
        etul_pnts.append(et_mjd[el])

    if ax_scale == 'lin':

      ax1.errorbar(t_pnts, pnts, xerr=et_pnts, yerr=e_pnts, fmt='o', mec='k', label='data')
      ax1.errorbar(tul_pnts, ul_pnts, xerr=[etul_pnts, etul_pnts], yerr=eul_pnts, uplims=True, fmt='bo', mec='k')

    if ax_scale == 'log' :
      ax2.errorbar(t_pnts, pnts, xerr=et_pnts, yerr=e_pnts, fmt='o', mec='k', label='data')
      ax2.errorbar(tul_pnts, ul_pnts, xerr=[etul_pnts, etul_pnts], yerr=eul_pnts, uplims=True, fmt='bo', mec='k')

  if ax_scale == 'lin':
    ax1.axis(axisLim) if axisLim != 'auto' else None
    ax1.grid()
    ax1.set_xlabel('t (MJD)')
    ax1.set_ylabel('dN/dE (ph/$cm^2$/s/MeV)')
    ax1.set_title('lightcurve') if title == 'none' else plt.title(title)
  #ax1.legend(loc=0)

  if ax_scale == 'log' :
    ax2.axis(axisLim) if axisLim != 'auto' else None
    ax2.grid()
    ax2.set_xlabel('t (MJD)')
    ax2.set_ylabel('dN/dE (ph/$cm^2$/s/MeV)')
    ax2.set_title('lightcurve') if title == 'none' else plt.title(title)
    #ax2.legend(loc=0)

  plt.subplots_adjust(hspace=0.5)

  # save fig ---!
  if ax_scale == 'lin' :
    extent = ax1.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(lightcurve.replace('.fits', '.png'), bbox_inches=extent.expanded(1.3, 1.3))

  if ax_scale == 'log' :
    extent = ax2.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(lightcurve.replace('.fits', '_log.png'), bbox_inches=extent.expanded(1.3, 1.3))

  # show fig ---!
  plt.show() if show == True else None
  plt.close()

  return

# WIP ---!
def showTSmap(tsmap, reg='none', col='black', suffix='none', title='none', show=True):
  '''
  :param:
  tsmap = TSV map fits file (str)
  show = if calling image show (bool)
  :return:
  '''

  # info data ---!
  hdulist = fits.open(tsmap)
  #print(hdulist.info())

  # get data ---!
  data = hdulist[0].data
  #print('\n shape:', data.shape, '\n type:', data.dtype.name)

  hdulist.close()

  # PLOT ---!
  ax = plt.subplot(111)
  plt.rc('text', usetex=True)

  # handle region ---!
  if reg != 'none':
    r = pyregion.open(reg).as_imagecoord(hdulist[0].header)
    for i in range(len(r)):
      r[i].attr[1]['color'] = col
      patch_list, text_list = r.get_mpl_patches_texts()
      for p in patch_list:
        ax.add_patch(p)
      for t in text_list:
        ax.add_artist(t)

  plt.imshow(data, cmap='bwr')
  #plt.xlabel('R.A. (deg)')
  #plt.ylabel('Dec (deg)')
  plt.xlabel('pix')
  plt.ylabel('pix')
  plt.title('TS map') if title == 'none' else plt.title(title)
  plt.colorbar().set_label('Significance TSV')

  # save fig ---!
  if suffix != 'none':
    plt.savefig(tsmap.replace('.fits', '_%s.png' % suffix))
  else:
    plt.savefig(tsmap.replace('.fits', '.png'))

  #show fig ---!
  plt.show() if show == True else None
  plt.close()

  return

# v02 ---!
def showButterfly_v02(diagram, spectrum, axisLim='auto', suffix='none', title='none', show=True):
  '''
  :param
  diagram = butterfly diagram txt file (str)
  spectrum = spectrum fits file (str)
  flux_pnts = list of fluxes in ph/cm^2/s/MeV (list of floats)
  fluxEn_pnts = list of energies in MeV
  show = if calling image show (bool)
  :return:
  '''
  data = np.loadtxt(diagram, delimiter=' ')

  hdulist = fits.open(spectrum)

  # get data ---!
  spc = hdulist[1].data

  # store data ---!
  en = spc.field(0)
  en_down = spc.field(1)
  en_up = spc.field(2)
  flux_pnt = spc.field(3)
  err_flux = spc.field(4)

  hdulist.close()

  # adjusting errors ---!
  en_err = (en_up - en_down) / 2

  # PLOT ---!
  plt.rc('text', usetex=True)

  # store data ---!
  energy = data[:, 0]
  intensity = data[:, 1]
  lerr = data[:, 2]
  uerr = data[:, 3]

  # from intensity (ph/cm^2/s/MeV) find flux in (erg/cm^2/s) ---!
  (flux, f_lerr, f_uerr) = ((intensity, lerr, uerr) * (energy ** 2)) / (6.4215 * 1e5)

  # PLOT ---!

  ax = plt.subplot(111, yscale='log', xscale='log')
  plt.rc('text', usetex=True)

  plt.plot(energy / 1e6, flux, 'b-', alpha=1, label='best fit')
  plt.fill_between(energy / 1e6, f_lerr, f_uerr, facecolor='blue', alpha=0.3)
  plt.scatter(en, flux_pnt, marker='+', c='r', label='spc', alpha=1)

  ax.axis(axisLim) if axisLim != 'auto' else None
  plt.xlabel('Energy (TeV)')
  plt.ylabel('Flux (erg/$cm^2$/s)')
  plt.title('butterfly diagram') if title == 'none' else plt.title(title)
  plt.legend(loc=0)
  plt.grid(True)

  # save fig ---!
  if suffix != 'none':
    plt.savefig(diagram.replace('.txt', '_%s.png' % suffix))
  else:
    plt.savefig(diagram.replace('.txt', '.png'))

  # show fig ---!
  plt.show() if show == True else None
  plt.close()

  return

# V01 ---!
def degradedIRF_3d(x, y, z, xlabel='x', ylabel='y', zlabel='z', title=None, c=['b'], zscale='linear',
                   fontsize=14, zlim=(0,1), alpha=[1], label=None, savefig=None, show=True) :

  fig = plt.figure(figsize=(12, 6))
  plt.rc('text', usetex=True)
  sns.set_style("whitegrid", {'axes.grid': False})
  ax = fig.add_subplot(111, projection='3d', zscale=zscale)

  curve = []
  for i in range(len(z)) :
    ax.plot_surface(x, y, z[i], alpha=alpha[i], color=c[i], label=label[i])
    curve.append(Rectangle((0, 0), 1, 1, fc=c[i], fill=True))

  ax.set_zlim(zlim)
  ax.set_xlabel(xlabel, fontsize=fontsize)
  ax.set_ylabel(ylabel, fontsize=fontsize)
  ax.set_zlabel(zlabel, fontsize=fontsize, labelpad=12)
  ax.set_title(title, fontsize=fontsize) if title != None else None
  ax.legend(curve, label, loc=0) if label != None else None

  plt.tight_layout()
  fig.savefig(savefig) if savefig != None else None
  plt.show() if show != False else None
  plt.close()

  return

# V01 ---!
def interp_ebl(x, y, savefig, kind='linear', xlabel='x', ylabel='y', title='title',
               label=['y', 'y2'], fontsize=12, show=True) :

  fig = plt.figure()
  plt.rc('text', usetex=True)
  sns.set()

  ax = plt.subplot(111, xscale='log', yscale='log')
  ax.plot(x[0], y[0], '.', label=label[0], c='g')
  ax.plot(x[1], y[1], 'o', c='k', markeredgecolor='k', markerfacecolor='none', label=label[1])
  ax.set_ylabel(ylabel, fontsize=fontsize)
  ax.set_xlabel(xlabel, fontsize=fontsize)
  ax.set_title(title, fontsize=fontsize)
  ax.legend(loc=0)

  plt.tight_layout()
  fig.savefig(savefig)
  plt.show() if show==True else None

  return

# SENSITIVITY ---!
def showSensitivity(x, y, savefig, xlabel='x', ylabel='y', label=['y'], title='none', fontsize=12, marker=['.'], show=True) :

  fig = plt.figure()
  plt.rc('text', usetex=True)
  sns.set()

  ax = plt.subplot(111, xscale='log', yscale='log')
  for i in range(len(y)) :
    ax.plot(x[i], y[i], marker=marker[i], label=label[i])
  ax.set_ylabel(ylabel, fontsize=fontsize)
  ax.set_xlabel(xlabel, fontsize=fontsize)
  ax.set_title(title, fontsize=fontsize) if title!='none' else None
  ax.legend(loc=0)

  plt.tight_layout()
  fig.savefig(savefig)
  plt.show() if show==True else None

  return
