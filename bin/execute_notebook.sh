#!/bin/bash

#########################################################
# USAGE:
#   execute_notebook.sh ${id} ${task.cpus} notebook.Rmd "-r xxx=yyy -r aaa=bbb"
#########################################################

export OPENBLAS_NUM_THREADS=$2 OMP_NUM_THREADS=$2  MKL_NUM_THREADS=$2 OMP_NUM_cpus=$2  MKL_NUM_cpus=$2 OPENBLAS_NUM_cpus=$2

jupytext -k python3 --to ipynb -o ${1}.raw.ipynb $3 \
  && papermill ${1}.raw.ipynb ${1}.ipynb $4 \
  && jupyter nbconvert --to html ${1}.ipynb
