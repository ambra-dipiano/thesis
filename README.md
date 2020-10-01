# thesis
Contains libs and scripts for my thesis project: detection of short gamma-ray bursts with CTA through real-time analysis. In order to run the scripts and use the defined classes, a config.xml file must be setup (see the "config.xml" example file) with at least one path for the "root" folder.

0. config.xml :: configuration file example
1. degradation.py :: script for calibration database degradation
2. modeule_plot.py :: functions for astronomical plots
3. module_statistics.py :: functions for statistical plots and basic analysis
4. pkg_blindsearch.py :: classes and functions for rta including setup, blind-search, data analysis and  model manipulation
5. rta-blindsrch.py :: script for rta statistical analysis of a given time frame
6. rta-followup.py :: script for rta follow-up of given observations
7. rta-wilks.py :: script for empty fields analysis in a given time frame
8. sensitivity.py :: script for sensitivity computation with a given IRF
9. templates.py :: script for plotting spectra and lightcurves of a given template
10. wilks.py :: script for post-processing analysis of empty fields

## Installation

Choose a virtual environment name and create it with the following command: 
```bash
conda env create --name \<envname\> --file=environments.yml
```
The above command will also install the required depedencies.

Test the ctools library with:
```bash
python -c 'import ctools; ctools.test()'
python -c 'import cscripts; cscripts.test()'
python -c 'import gammalib; gammalib.test()'
```


### cfitsio not found bug fix
ctools 1.7.1 will look for cfitsio.so.5 but the current version of the anaconda package contains cfitsio.so.8 . 
If the cfitsio library is not found, issue the following command:
```bash
mv \<environment path\>/lib/libcfitsio.so.8 \<environment path\>/lib/libcfitsio.so.5
```