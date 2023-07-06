#!/usr/bin/bash


cd /share/pub/likai/Star_protocols/scripts/03
for file in 03_11_HardFiltering_chr*.sh
do
        echo $file
        sbatch --mem=20G $file
done
