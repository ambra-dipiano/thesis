>> DIRECTORIES:
!! run0406 ---> run0406_ID000126 L.Nava GRB template from GammaCatalogV1.0

>> FILES:
!! script_RTAdetection_v05.py ---> simulation of a chunk for all texp=1,5,10,100s, from template model + detection chain, results written in csv file
!! run0406_task001.sh - test script and job submission (10 ph-lists simulations only)
!! script_RTAdetection_v06.py ---> simulation of a chunk and fixed texp, from template model + detection chain, results written in csv file
!! module_xml.py ---> collection of useful python function, among which the ones used to perform dection, max likelihood etc
!! run0406_task001.sh ---> testing script for job submission and script

>> TAR.GZ:
!! run0406_10x10s.tar.gz ---> chunk 1 of 10 sim for texp=10s (trial to solve lxml bug) [FIXED]


>> TO-DO:
1) test first chunk of simulations & script: DONE
2) improve csv writing as to avoid over-writing different chunks (try ... except ...): DONE
3) test global count for chunks and simulations id: DONE
4) script to run parallel chunks (100 trials each for texp=1,5,10,100s) but not ctobssim: DONE
5) filter Nsrc e TSV etc, so that if TS < 9 skips: DONE
6) each chunk for all texps: DONE
