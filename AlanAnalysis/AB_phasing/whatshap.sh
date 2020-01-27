#!/bin/bash

ijob -c1 -p standard -A berglandlab
module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6

### install whatsapp
  #conda create -y -n whatshapp-env

  #source activate whatshapp-env
  #conda config --add channels bioconda

  #conda install whatshap nomkl

### load env
  source activate whatshapp-env
  whatshap_dir=/home/aob2x/.conda/pkgs/whatshap-0.18-py36h6bb024c_0/bin

python ${whatshap_dir}/whatshap
