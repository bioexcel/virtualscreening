export PYMDSETUP=$HOME/projects/pymdsetup
export PATH=$PYMDSETUP/gromacs_wrapper:$PYMDSETUP/gnuplot_wrapper:$PYMDSETUP/gromacs_wrapper:$PYMDSETUP/scwrl_wrapper:$PATH

cwl-runner scwrl.cwl test/scwrl_conf.yml
cwl-runner pdb2gmx.cwl test/pdb2gmx_conf.yml
cwl-runner editconf.cwl test/editconf_conf.yml
cwl-runner solvate.cwl test/solvate_conf.yml
cwl-runner grompp.cwl test/gppions_conf.yml
cwl-runner genion.cwl test/genion_conf.yml
cwl-runner mdrun.cwl test/mdmin_conf.yml
cwl-runner rms.cwl test/rms_conf.yml
cwl-runner gnuplot.cwl test/gnuplot_conf.yml
rm sol.gro sol_top.zip p2g.gro p2g_top.zip ec.gro gio.gro gio_top.zip gppions.tpr mdmin.gro mdmin.trr gplot.png mutated.pdb rms.xvg 2> /dev/null
