export PYMDSETUP=$HOME/projects/pymdsetup
export PATH=$PYMDSETUP/gromacs_wrapper:$PYMDSETUP/gnuplot_wrapper:$PYMDSETUP/gromacs_wrapper:$PYMDSETUP/scwrl_wrapper:$PATH

cwl-runner mutations.cwl mutations_conf.yml
