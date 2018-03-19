#!/bin/bash

enqueue_compss \
  --job_dependency=$1 \
  --exec_time=$2 \
  --num_nodes=$3 \
  --max_tasks_per_node=1 \
  --qos=bsc_ls \
  --worker_working_dir=gpfs \
  --network=infiniband \
  --lang=python \
  --pythonpath=/gpfs/home/bsc23/bsc23210/pymdsetup/:/gpfs/home/bsc23/bsc23210/ \
  --master_working_dir=/gpfs/scratch/bsc23/bsc23210/ \
  --worker_working_dir=/gpfs/scratch/bsc23/bsc23210/ \
  --tracing=$4 \
  --graph=$5 \
  --log_level=debug \
/gpfs/home/bsc23/bsc23210/pymdsetup/workflows/gromacs_full_pycompss.py $6 $7 $3 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17} ${18} ${19} ${20}

#./run_pycompss_mn.sh None 30 3 false false workflows/conf_2mut_nt0.yaml mare_nostrum
