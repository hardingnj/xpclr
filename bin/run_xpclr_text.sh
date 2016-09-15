#!/bin/bash

/home/njh/pyenv/python34/bin/python ~/git/xpclr/compute_xpclr.py \
  --map /home/njh/mapfile.snp \
  --popA /home/njh/genotype1.geno \
  --popB /home/njh/genotype2.geno \
  -C 3L \
  -O /tmp/BFS_vs_GWA_3L_text.txt
