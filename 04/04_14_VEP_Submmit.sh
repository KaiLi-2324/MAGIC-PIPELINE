#!/usr/bin/bash


cd /share/pub/likai/Star_protocols/scripts/04
for file in $(ls 04_14_VEP_Chr*)
do
        echo $file
        sbatch --mem=20G $file
done
