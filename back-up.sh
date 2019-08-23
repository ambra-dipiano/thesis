#!/bin/bash

<<<<<<< HEAD
rm run0406*task*.sh
rm job_run.sh
cp *py *sh bkup_scripts/.
cp $HOME/*py $HOME/*sh bkup_scripts/.
cp -r bkup_scripts $HOME/.
# ------------------------------------------------------------------------ TESTS !!!
#tar cvfz run0406_bkup.tar.gz ./run0406
#mv run0406_bkup.tar.gz $HOME/.
# ------------------------------------------------------------------------
#tar cvfz run0406_1000x_flux.tar.gz ./run0406 
#mv run0406_1000x_flux.tar.gz $HOME/.
# ------------------------------------------------------------------------
tar cvfz run0406_test03_binsky.tar.gz ./run0406
mv run0406_test03_binsky.tar.gz $HOME/.
# ------------------------------------------------------------------------
#tar cvfz crab_2000s_test05.tar.gz ./crab 
#mv crab_2000s_test05.tar.gz $HOME/.
# ------------------------------------------------------------------------ EBL !!!
#tar cvfz template_ebl_test01.tar.gz ./template_ebl_test
#mv template_ebl_test01.tar.gz $HOME/.
# ------------------------------------------------------------------------ WILKS !!!
#tar cvfz run0406_bkg_02.tar.gz ./run0406_bkg
#mv run0406_bkg_02.tar.gz $HOME/.
# ------------------------------------------------------------------------ LIGHTCURVES !!!
#tar cvfz crab_lc_v03.tar.gz ./crab_lc 
#mv crab_lc_v03.tar.gz $HOME/.
# --------------------
#tar cvfz run0406_lc_v03.tar.gz ./run0406_lc 
#mv run0406_lc_v03.tar.gz $HOME/.
# --------------------
=======
scp -r pianoambra@morgana.iasfbo.inaf.it:~/bkup_scripts .
cp bkup_scripts/module*py .
#scp pianoambra@morgana.iasfbo.inaf.it:~/run0406_test03_binsky.tar.gz . 
#scp pianoambra@morgana.iasfbo.inaf.it:~/template_ebl_test01.tar.gz .
#scp pianoambra@morgana.iasfbo.inaf.it:~/run0406_bkup.tar.gz . 
#scp pianoambra@morgana.iasfbo.inaf.it:~/run0406_1000x_flux.tar.gz . 
#scp pianoambra@morgana.iasfbo.inaf.it:~/crab_2000s_test05.tar.gz .
#scp pianoambra@morgana.iasfbo.inaf.it:~/run0406_bkg_02.tar.gz .
#scp pianoambra@morgana.iasfbo.inaf.it:~/crab_lc_v02.tar.gz .
#scp pianoambra@morgana.iasfbo.inaf.it:~/run0406_lc_v02.tar.gz .

>>>>>>> 620d34c4f29370e844e6d3556dc5c307ce76c438

