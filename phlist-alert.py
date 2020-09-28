# IMPORTS ---!
from pkg_blindsearch import *
import os

# import time

# --------------------------------- SETUP --------------------------------- !!!

# compact initialisation ---!
trials = 1  # trials
count = 0  # starting count
# cpus ---!
nthreads = 2
os.environ['OPENBLAS_NUM_THREADS'] = str(nthreads)
os.environ['MKL_NUM_THREADS'] = str(nthreads)

# ctools/cscripts parameters ---!
run_id = 'run0367'
merger_id = 'ID000227'
caldb = 'prod3b-v2'
irf = 'South_z40_0.5h'
talert = 120  # alert latency (s)
tslew = 30  # slewing time (s)
tdelay = talert + tslew  # total time wrt burst onset (s)
ttotal = 1e2  # total observation time (s)
run_duration = 1200  # observation run time (s) ---!
emin = 0.03  # simulation minimum energy (TeV)
emax = 150.0  # simulation maximum energy (TeV)
roi = 5  # region of interest for simulation and selection (deg)

# conditions control ---!
if_ebl = True  # uses the EBL absorbed template
add_ebl = True  # add ebl to template
extract_spec = True  # generates spectral tables and obs definition models
irf_degrade = False  # use degraded irf
use_run = False

# path configuration ---!
cfg = xmlConfig()
p = ConfigureXml(cfg)
p.setRunDir('%s_%s/' % (run_id, merger_id))
# files ---!
source_name = '%s_%s' % (run_id, merger_id)
ebl_table = p.getRootDir() + '/ebl_tables/gilmore_tau_fiducial.csv'
merger_map = p.getMergersDir() + '%s_Merger%s_skymap.fits' % (run_id, merger_id)
template = p.getTemplatesDir() + '%s.fits' % source_name
model_pl = p.getModelsDir() + '%s.xml' % source_name
tcsv = p.getDataDir() + 'time_slices.csv'

# check merger map existance
if not os.path.isfile(template):
    raise ValueError('No template is available for %s_%s' % (run_id, merger_id))
if not os.path.isfile(merger_map + '.gz'):
    raise ValueError('No merger map is available for %s_%s' % (run_id, merger_id))

# pointing with off-axis equal to max prob GW ---!
true_coord, pointing, offmax = getPointing(fits_file=template, merger_map=merger_map)
print('coords true:', true_coord, 'point', pointing, 'off', offmax)

# check source xml model template and create if not existing ---!
if not os.path.isfile(model_pl):
    print('Creating template model xml')
    template_pl = p.getModelsDir() + 'grb_file_model.xml'
    os.system('cp %s %s' % (template_pl, model_pl))
    model_xml = ManageXml(xml=model_pl)
    model_xml.setModelParameters(source=source_name, parameters=('RA', 'DEC'), values=true_coord)
    del model_xml

# recap and dof ---!
dof, m2, m1 = getDof()
print('!!! *** !!! dof = ', m2, ' - ', m1, ' = ', dof)
print('!!! *** !!! EBL ABSORPTION:', if_ebl)
print('!!! *** !!! IRF DEGRADATION:', irf_degrade)
print('!!! *** !!! nominal caldb:', caldb)
print('!!! *** !!! irf:', irf)
print('!!! *** !!! sim energy range: [', emin, ', ', emax, '] (TeV)')
print('!!! *** !!! roi:', roi, ' (deg)')
print('!!! *** !!! pointing:', pointing, ' (deg)')
print('!!! *** !!! delay time:', tdelay, ' s')
print('!!! *** !!! total observation time:', ttotal - tdelay, ' s')

# --------------------------------- INITIALIZE --------------------------------- !!!

# setup trials obj ---!
simulation = Analysis()
simulation.nthreads = nthreads
simulation.pointing = pointing
simulation.roi = roi
simulation.e = [emin, emax]
simulation.tmax = ttotal
simulation.model = model_pl
# degrade IRF if required ---!
simulation.caldb = caldb
simulation.irf = irf
simulation.template = template
# add EBL to template ---!
if add_ebl:
    print('Computing EBL absorption')
    simulation.table = ebl_table  # fiducial ---!
    simulation.zfetch = True
    simulation.if_ebl = False
    simulation.fitsEbl(template.replace('.fits', '_ebl.fits'), ext_name='EBL-ABS. SPECTRA')
    simulation.if_ebl = if_ebl
    simulation.template = template.replace('.fits', '_ebl.fits')
print('!!! check ---- template:', simulation.template)
# load template ---!
simulation.extract_spec = extract_spec
tbin_stop = simulation.loadTemplate(source_name=source_name, return_bin=True, data_path=p.getDataDir())
print('!!! check ---- tbin_stop:', tbin_stop)
print('!!! check ---- caldb:', simulation.caldb)

# --------------------------------- removed 1° LOOP :: trials  --------------------------------- !!!

count += 1
simulation.seed = count
# attach ID to fileroot ---!
if irf_degrade:
    prefix = 'deg'
else:
    prefix = 'ebl'
f = '%s%06d' % (prefix, count)
print('!!! check ---- grb root name:', f)

event_bins = []
# --------------------------------- SIMULATION GRB --------------------------------- !!!
simulation.table = tcsv
tgrid = simulation.getTimeSlices(GTI=(tdelay, ttotal))  # methods which returns time slice edges ---!
# simulate ---!
for i in range(tbin_stop):
    simulation.t = [tgrid[i], tgrid[i + 1]]
    if if_ebl:
        simulation.model = p.getDataDir() + '%s_ebl_tbin%02d.xml' % (source_name, i)
    else:
        simulation.model = p.getDataDir() + '%s_tbin%02d.xml' % (source_name, i)
    event = p.getObsDir() + "%s_tbin%02d.fits" % (source_name, i)
    event_bins.append(event)
    simulation.output = event
    simulation.eventSim()

if not use_run:
    # --------------------------------- APPEND EVENTS IN SINGLE PH-LIST --------------------------------- !!!
    event_list = p.getObsDir() + f + '.fits'
    simulation.input = event_bins
    simulation.output = event_list
    simulation.appendEventsSinglePhList(GTI=[tdelay, ttotal])
    phlist = event_list
    print('phlist name', phlist)
else:
    # --------------------------------- APPEND EVENTS IN MULTIPLE PH-LIST --------------------------------- !!!
    event_list = p.getObsDir() + f + '.fits'
    simulation.input = event_bins
    simulation.output = event_list
    num_max, phlist = simulation.appendEventsMultiPhList(max_length=run_duration, last=ttotal)
    print('runs (', num_max, ') name', phlist)

del simulation
print('remove ' + p.getObsDir() + '*tbin*')
os.system('rm ' + p.getObsDir() + '*tbin*')
