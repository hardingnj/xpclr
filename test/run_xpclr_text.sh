#!/bin/bash

python ../bin/compute_xpclr.py \
  --map ../example/mapfile.snp \
  --popA ../example/genotype1.geno \
  --popB ../example/genotype2.geno \
  -C 1 \
  -O ../example/output.txt
