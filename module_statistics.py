# ============================================
# !!! MODULE FOR STATISTICS AND HISTOGRAMS !!!
# ============================================

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.patches import Rectangle
from scipy import stats
from scipy.stats import rayleigh, norm, chi2
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Circle

extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
extra2 = Line2D([0], [0], ls='-.', color='k', lw='1')


# HIST 1D GAUSSIAN DISTRIBUTION ---!
def hist1d_gauss(x, mean, loc=0, threshold=1, nbin=20, width=None, fontsize=12, color='b', alpha=0.5, title='gaussian fit',
                 ax_thresh=0.2, xlabel='x', ylabel='y',filename='hist1d_gauss.png') :

    if width == None :
        width = threshold/nbin
    if nbin == None :
        nbin = int(threshold/width)
    if nbin == None and width == None :
        print('Error: set either nbin or width')

    fig = plt.figure()
    plt.rc('text', usetex=True)
    sns.set()

    ax = plt.subplot(111)
    sns.distplot(x, bins=nbin, kde=False,
                 fit=norm, norm_hist=True,
                 color=color, hist_kws={'alpha':alpha, 'range':[loc-threshold, loc+threshold]})
    plt.plot([],[], color='none', label='wbin=%.2fdeg' %width)
    plt.axvline(loc, c='r', ls='-.', label='true $\\approx$ %.3fdeg' %loc)
    plt.axvline(mean, c='k', ls='--', label='mean $\\approx$ %.3fdeg' %mean)
    plt.title(title, fontsize=fontsize)
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.legend(fontsize=fontsize)
    plt.xlim([loc-ax_thresh, loc+ax_thresh])

    fig.savefig(filename)

    return fig, ax


# HIST 1D RAYLEIGH DISTRIBUTION ---!
def hist1d_rayleigh(x, mean, rayleigh_prms={'loc':0, 'scale':1}, threshold=1, nbin=None, width=None, fontsize=12, color='b', alpha=0.5, title='rayleigh fit',
                    ax_thresh=0.2, xlabel='x', ylabel='y', filename='hist1d_rayleigh.png') :

#    MLE = np.exp(-0.5)/rayleigh_prms['scale'] if rayleigh_prms != None else None

    if width == None :
        width = threshold/nbin
    if nbin == None :
        nbin = int(threshold/width)
    if nbin == None and width == None :
        print('Error: set either nbin or width')

    fig = plt.figure()
    plt.rc('text', usetex=True)
    sns.set()

    ax = plt.subplot(111)
    sns.distplot(x, bins=nbin, kde=False,
                 fit=rayleigh, norm_hist=True,
                 color=color, hist_kws={'alpha':alpha, 'range':[0.0, threshold]})
    plt.plot([],[], color='none', label='wbin=%.2fdeg' %width)
    plt.axvline(rayleigh_prms['scale'], c='r', ls='-.', label='mode $\\approx$ %.3fdeg' %rayleigh_prms['scale']) if rayleigh_prms['scale'] != None else None
#    plt.axvline(MLE, c='m', ls='-.', label='MLE $\\approx$ %.3fdeg' %MLE) if rayleigh_prms != None else None
    plt.axvline(mean, c='k', ls='--', label='mean $\\approx$ %.3fdeg' %mean)
    plt.title(title, fontsize=fontsize)
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.legend(fontsize=fontsize)
    plt.xlim([rayleigh_prms['loc'], rayleigh_prms['loc']+ax_thresh]) if rayleigh_prms['loc'] != None else None

    fig.savefig(filename)

    return fig, ax


# RAYLEIGH CDF WITH CONFIDENCE INTERVAL ---!
def rayleigh_cdf(x, loc=0, scale=1, if_CI=True, probs=(0.6827, 0.9545, 0.9973, 0.99994), xlabel='x', title='x$\\sim$ RA($\\gamma$) CDF',
                    colors=('k', 'r', 'orange', 'm'), fontsize=12, filename='theo_rayleigh_cdf.png') :

    fig = plt.figure()
    plt.rc('text', usetex=True)
    sns.set()

    ax = plt.subplot(111)
    ax.plot(np.sort(x), stats.rayleigh.cdf(np.sort(x), loc=loc, scale=scale), ls='-', label='cdf')
    ax.axvline(scale, c='maroon', label='$\gamma$')
    ax.axvline(np.std(x), c='maroon', ls=':', label='1 std =%.2f' %(np.std(x)))

    if if_CI is True :
        x_critical = []
        for i in range(len(probs)) :
            x_critical.append(stats.rayleigh.ppf(q=probs[i], loc=loc, scale=scale))
            ax.axvline(x_critical[i], c=colors[i], ls='-.',
                         label='x=%.2f, %.2f' %(x_critical[i],probs[i]*100)+'%')

    plt.ylabel('1-$\\alpha$', rotation=90, fontsize=fontsize)
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.title(title, fontsize=fontsize)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    plt.legend(loc=0)

    fig.savefig(filename)

    return fig, ax


# RAYLEIGH PDF WITH CONFIDENCE INTERVAL ---!
def rayleigh_pdf(x, loc=0, scale=1, if_CI=True, probs=(0.6827, 0.9545, 0.9973, 0.99994), xlabel='x', title='x$\\sim$ RA($\\gamma$) CDF',
                    colors=('k', 'r', 'orange', 'm'), fontsize=12, filename='theo_rayleigh_cdf.png') :

    fig = plt.figure()
    plt.rc('text', usetex=True)
    sns.set()

    ax = plt.subplot(111)
    ax.plot(np.sort(x), stats.rayleigh.pdf(np.sort(x), loc=loc, scale=scale), ls='-', label='cdf')
    ax.axvline(scale, c='maroon', label='$\gamma$')
    ax.axvline(np.std(x), c='maroon', ls=':', label='1 std =%.2f' %(np.std(x)))

    if if_CI is True :
        x_critical = []
        for i in range(len(probs)) :
            x_critical.append(stats.rayleigh.ppf(q=probs[i], loc=loc, scale=scale))
            ax.axvline(x_critical[i], c=colors[i], ls='-.',
                         label='x=%.2f, %.2f' %(x_critical[i],probs[i]*100)+'%')

    plt.ylabel('counts density', rotation=90, fontsize=fontsize)
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.title(title, fontsize=fontsize)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    plt.legend(loc=0)

    fig.savefig(filename)

    return fig, ax


# 2D HISTOGRAM WITH RAYLEIGH CONFIDENCE INTERVAL ---!
def hist2d_rayleigh_CI(x, y, nbin=None, width=None, rayleigh_prms={'loc':0, 'scale':1}, xcentre=0, ycentre=0, threshold=1, probs=(0.6827, 0.9545, 0.9973, 0.99994), colors=('k', 'r', 'orange', 'm'), ax_thresh=0.2, xlabel='x', ylabel='y', title='confidence intervals from theoretical distribution', fontsize=12, filename='hist2d_CIrayleigh.png') :

    xmean = np.mean(x)
    ymean = np.mean(y)

    if width is None :
        width = threshold/nbin
    if nbin is None :
        nbin = int(threshold/width)
    if nbin is None and width is None :
        print('Error: set either nbin or width')

    fig = plt.figure()
    plt.rc('text', usetex=True)
    sns.set()

    ax = plt.subplot(111)
    h = ax.hist2d(x, y, bins=nbin, cmap='jet',
                  range=[[xcentre - threshold, xcentre + threshold], [ycentre - threshold, ycentre + threshold]])
    plt.scatter(xcentre, ycentre, c='w', marker='*', s=1e2)
    plt.plot([], [], c='none', label='Reyleigh')
    for i in range(len(probs)):
        plt.plot([], [], c=colors[i], label='%.2f' % (probs[i] * 100) + '\%')
        r = stats.rayleigh.ppf(q=probs[i], loc=rayleigh_prms['loc'], scale=rayleigh_prms['scale'])
#        q = rayleigh['scale'] * np.sqrt(-2 * np.log(probs[i]))
#        r = stats.rayleigh.ppf(q=q, loc=rayleigh['loc'], scale=rayleigh['scale'])
        cir = Circle(xy=(xmean, ymean),
                     radius=r,
                     color=colors[i], lw=2)
        cir.set_facecolor('none')
        ax.add_artist(cir)

    plt.colorbar(h[3], ax=ax).set_label('counts')
    plt.axis([xcentre - ax_thresh, xcentre + ax_thresh, ycentre - ax_thresh, ycentre + ax_thresh], 'equal') if ax_thresh != None else None
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.title(title, fontsize=fontsize)
    plt.legend(ncol=3)

    fig.savefig(filename)

    return fig, ax


# COVARIANCE EIGENVALUES ---!
def eigsorted(cov):
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]

    return vals[order], vecs[:, order]


# 2D HISTOGRAM WITH GAUSSIAN COVARIANCE CONFIDENCE INTERVAL ---!
def hist2d_gauss_CI(x, y, nbin=None, width=None, xcentre=0, ycentre=0, threshold=1, nstd=(1, 2, 3, 5),
                    colors=('k', 'r', 'orange', 'm'), ax_thresh=0.2, xlabel='x', ylabel='y',
                    title='confidence intervals from theoretical distribution', fontsize=12, filename='hist2d_CIgauss.png') :

    xmean = np.mean(x)
    ymean = np.mean(y)

    if width is None :
        width = threshold/nbin
    if nbin is None :
        nbin = int(threshold/width)
    if nbin is None and width is None :
        print('Error: set either nbin or width')

    fig = plt.figure()
    plt.rc('text', usetex=True)
    sns.set()

    ax = plt.subplot(111)
    h = ax.hist2d(x, y, bins=nbin, cmap='jet',
                   range=[[xcentre - threshold, xcentre + threshold], [ycentre - threshold, ycentre + threshold]])
    plt.scatter(xcentre, ycentre, c='w', marker='*', s=1e2)
    plt.plot([], [], c='none', label='gauss')
    for i in range(len(nstd)):
        plt.plot([], [], c=colors[i], label='%d $\sigma$' % (nstd[i]))
        cov = np.cov(x, y)
        vals, vecs = eigsorted(cov)
        theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
        w, v = 2 * nstd[i] * np.sqrt(vals)
        ell = Ellipse(xy=(xmean, ymean),
                      width=w, height=v,
                      angle=theta, color=colors[i], lw=2)
        ell.set_facecolor('none')
        ax.add_artist(ell)

    plt.colorbar(h[3], ax=ax).set_label('counts')
    plt.axis([xcentre - ax_thresh, xcentre + ax_thresh, ycentre - ax_thresh, ycentre + ax_thresh], 'equal') if ax_thresh != None else None
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.title(title, fontsize=fontsize)
    plt.legend(ncol=3)

    fig.savefig(filename)

    return fig, ax


# 2D HISTOGRAM MAP ---!
def hist2d_map(x, y, nbin=None, width=None, xcentre=0, ycentre=0, threshold=1, ax_thresh=0.2, xlabel='x', ylabel='y',
              title='probability map', fontsize=12, filename='hist2d_map.png', if_CI=None, rayleigh={'loc':0, 'scale':1},
              nstd=(1, 2, 3, 5), colors=('k', 'r', 'orange', 'm'), probs=(0.6827, 0.9545, 0.9973, 0.99994), smooth=True) :

    if width is None :
        width = threshold/nbin
    if nbin is None :
        nbin = int(threshold/width)
    if nbin is None and width is None:
        print('Error: set either nbin or width')

    fig = plt.figure()
    plt.rc('text', usetex=True)
    sns.set()

    ax = plt.subplot(111)
    h = ax.hist2d(x, y, bins=nbin, cmap='jet', density=True,
                  range=[[xcentre - threshold, xcentre + threshold], [ycentre - threshold, ycentre + threshold]])
    if smooth:
        plt.clf()
        plt.close()
        plt.imshow(h[0], origin="lower", interpolation="gaussian")
    plt.scatter(xcentre, ycentre, c='w', marker='*', s=1e2)

    if if_CI is None :
        pass

    elif if_CI.lower() is 'rayleigh' :
        xmean = np.mean(x)
        ymean = np.mean(y)
        plt.plot([], [], c='none', label='Reyleigh')
        for i in range(len(probs)):
            plt.plot([], [], c=colors[i], label='%.2f' % (probs[i] * 100) + '%')
            r = stats.rayleigh.ppf(q=probs[i], loc=rayleigh['loc'], scale=rayleigh['scale'])
            cir = Circle(xy=(xmean, ymean),
                         radius=r,
                         color=colors[i], lw=2)
            cir.set_facecolor('none')
            ax.add_artist(cir)

    elif if_CI.lower() is 'gauss' or if_CI.lower is 'covariance' or if_CI.lower is 'cov' :
        xmean = np.mean(x)
        ymean = np.mean(y)
        plt.plot([], [], c='none', label='gauss')
        for i in range(len(nstd)):
            plt.plot([], [], c=colors[i], label='%.2f' % (nstd[i] * 100) + '%')
            cov = np.cov(x, y)
            vals, vecs = eigsorted(cov)
            theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
            w, v = 2 * nstd[i] * np.sqrt(vals)
            ell = Ellipse(xy=(xmean, ymean),
                          width=w, height=v,
                          angle=theta, color=colors[i], lw=2)
            ell.set_facecolor('none')
            ax.add_artist(ell)
    else :
        print('Error: if_CI parameter value not understood')

    plt.colorbar(h[3], ax=ax)
    plt.axis([xcentre - ax_thresh, xcentre + ax_thresh, ycentre - ax_thresh, ycentre + ax_thresh], 'equal')
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.title(title, fontsize=fontsize)

    fig.savefig(filename)

    return fig, ax


# WILKS THEOREM DIST FOR EMPTY FIELDS ---!
def ts_wilks(x, trials, df=1, nbin=None, width=None, ylim=None, xlim=None, xlabel='TS',
             fontsize=12, figsize=(4,6),
             title='TS distribution (empty fields)', filename='wilks_preTrials.png') :

    if width is None :
        width = (max(x)-min(x))/nbin
    if nbin is None :
        nbin = int((max(x)-min(x))/width)
    if nbin is None and width is None:
        print('Error: set either nbin or width')

    fig = plt.figure(figsize=figsize)
    plt.rc('text', usetex=True)
    sns.set()

    ax = plt.subplot(111, yscale='log')
    h, edges = np.histogram(x, bins=int(nbin), density=False, range=(0., max(x)))
    h_norm = h / trials
    cbin = (edges[1:] + edges[:-1]) / 2
    yerr = np.sqrt(h) / trials
    xerr = (edges[:-1] - edges[1:]) / 2

    x2 = []
    for i in range(len(x)):
        if x[i] > edges[1] :
            x2.append(x[i])

    plt.errorbar(cbin, h_norm, fmt='k+', yerr=yerr, xerr=xerr, markersize=5, label='data')
    plt.plot(x2, stats.chi2.pdf(x2, df=df) / 2, lw=1, ls='--', label='$\\chi^2$/2')

    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel('normalized counts', fontsize=fontsize)
    plt.title(title, fontsize=fontsize)
    plt.legend(loc=0, fontsize=fontsize)
    plt.xlim(xlim) if xlim is not None else None
    plt.ylim(ylim) if ylim is not None else None

    fig.savefig(filename)

    return fig, ax


# WILKS THEOREM P-VALUES FOR EMPTY FIELDS ---!
def p_values(x, trials, df=1, nbin=None, width=None, ylim=None, xlim=None, xlabel='h',
             fontsize=12, figsize=(4,6),
             title='p-value (empty fields)', filename='pvalue_preTrials.png') :

    if width is None :
        width = (max(x)-min(x))/nbin
    if nbin is None :
        nbin = int((max(x)-min(x))/width)
    if nbin is None and width is None:
        print('Error: set either nbin or width')

    fig = plt.figure(figsize=figsize)
    plt.rc('text', usetex=True)
    sns.set()

    ax = plt.subplot(111, yscale='log')

    h, edges = np.histogram(x, bins=int(nbin), density=False, range=(0., max(x)))
    yerr = np.sqrt(h) / trials
    h = h/trials
    cbin = (edges[1:] + edges[:-1]) / 2
    p = (1 - stats.chi2.cdf(cbin, df=df))
    xerr = (edges[:-1] - edges[1:]) / 2

    x2 = []
    for i in range(len(x)):
        if x[i] > edges[1] :
            x2.append(x[i])

    plt.errorbar(cbin[0], h[0], yerr=yerr[0], xerr=xerr[0], fmt='k+', markersize=5)
    plt.errorbar(cbin[1:], p[1:], yerr=yerr[1:], xerr=xerr[1:], fmt='k+', markersize=5, label='p-values')
    plt.plot(x2, (1 - stats.chi2.cdf(x2, df=1)), lw=1, ls='-.', c='maroon', label='$P$')

    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel('probability', fontsize=fontsize)
    plt.title(title, fontsize=fontsize)
    plt.legend(loc=0, fontsize=fontsize)
    plt.xlim(xlim) if xlim is not None else None
    plt.ylim(ylim) if ylim is not None else None

    fig.savefig(filename)

    return fig, ax
