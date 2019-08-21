from module_xml import *
from module_analysis import *
from module_plot import *

workdir = '/home/ambra/Desktop/cluster-morgana/run0406_plots/run0406/' 
path = workdir + 'run0406_ID000126/'
datapath = path + 'data/'
simpath = path + 'sim/'
selectpath = path + 'selected_sim/'
detpath = path + 'detection_all/'
template = workdir + 'run0406_ID000126.fits'
model = workdir + 'run0406_ID000126.xml'

count = 1
tmax= [150]
fileroot = 'run0406_'
f = fileroot + 'sim%06d' % (count)

t, tbin_stop = load_template(template, model, tmax, extract_spec=True, pathout=datapath) 

# ctools/cscripts parameters ---!
caldb = 'prod3b'
irf = 'South_z40_average_100s'

sigma = 5  # detection acceptance (Gaussian)
texp = [1, 5, 10, 100]  # exposure times (s) 
texp.sort()
tint = len(texp) 
tmin = 30  # slewing time (s)
tmax = []
for i in range(tint):
    tmax.append(tmin + texp[i])

elow = 0.03  # simulation minimum energy (TeV)
ehigh = 0.5  # simulation maximum energy (TeV)
emin = 0.03  # selection minimum energy (TeV)
emax = 0.5  # selection maximum energy (TeV)
roi = 5  # region of interest for simulation and selection (deg)
wbin = 0.02  # skymap bin width (deg)
nbin = int(roi/wbin)  # skymap x,y axes (pixel)
confidence = (0.68, 0.95, 0.9973)  # confidence interval for asymmetrical errors (%)

# pointing with off-axis equal to max prob GW ---!
offmax = (-1.475, -1.37)  # (deg)
trueRa = 33.057  # (deg)
trueDec = -51.841  # (deg)
pointRA = trueRa + offmax[0] # (deg)
pointDEC = trueDec + offmax[1] # (deg)

# ======================
# !!! GRB SIMULATION !!!
# ======================

event_bins = []

# simulate ---!
for i in range(tbin_stop):
    model = datapath + 'run0406_ID000126_tbin%d.xml' % i
    event = simpath + f + "_tbin%02d.fits" % i
    event_bins.append(event)
    if not os.path.isfile(event) :
        simulate_event(model=model, event=event, t=[t[i], t[1+i]], tbin_stop=tbin_stop, e=[elow, ehigh], caldb=caldb, irf=irf, pointing=[pointRA, pointDEC], seed=count) 

# observation list ---!
eventList = simpath + 'obs_%s.xml' % f
if not os.path.isfile(eventList) :
    observation_list(event=event_bins, eventList=eventList, obsname=f, tbin_stop=tbin_stop) 

print('!!! check --- list obs:', eventList)

# ============================
# !!! EVENT LIST SELECTION !!!
# ============================

selectedEvents = []

for i in range(tint) :
    selectedEvents.append(eventList.replace(simpath, selectpath).replace('obs_', 'texp%ds_' % texp[i]))
    prefix = selectpath + 'texp%ds_' % texp[i]
    if not os.path.isfile(selectedEvents[i]) :
        select_event(eventList=eventList, event_selected=selectedEvents[i], 
                     prefix=prefix, t=[tmin, tmax[i]], e=[emin, emax])   

print('!!! check --- selection: ', selectedEvents)

#os.system('python ./cta_scripts/cta_obs.py '+eventList)

# ========================
# !!! SELECTION SKYMAP !!!
# ========================

skymapName = []

for i in range(tint) :
    skymapName.append(selectedEvents[i].replace(selectpath, detpath).replace('.xml', '_skymap.fits'))
    if not os.path.isfile(skymapName[i]) :
        skymap_event(event_selected=selectedEvents[i], sky=skymapName[i], e=[emin, emax], 
                     caldb=caldb, irf=irf, wbin=wbin, nbin=nbin)

print('!!! check --- skymaps: ', skymapName)


for i in range(tint) :
    showSkymap(skymapName[i], title='texp=%ds' %texp[i], show=False)

# ==============================
# !!! DETECTION AND MODELING !!!
# ==============================

detXml = []
detReg = []
pos = []

for i in range(tint) :
    det, reg, coord = srcDetection_spcModeling(skymapName[i], sigma=sigma, maxSrc=10)
    detXml.append(det)
    detReg.append(reg)
    pos.append(coord)


    print('\n\n==========\n\n!!! check --- detection.............', texp[i],'s done')
    print('\n\n!!! coords:', pos[i], '\n\n ==========\n\n', detReg[i])

for i in range(tint) :
    showSkymap(skymapName[i], reg=detReg[i], col='black', suffix='det%dsgm_texp%ds' %(sigma,texp[i]), title='detection sgm=%d, texp=%ds' %(sigma,texp[i]), show=False)
