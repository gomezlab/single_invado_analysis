#!/usr/bin/env bash
echo "Syncing Focal Adhesion Results to $1";
rsync --progress -a ../../results/focal_adhesions/ $1:~/Documents/Projects/focal_adhesions/results/focal_adhesions/;

echo "Syncing S178A Results to $1";
rsync --progress -a ../../results/S178A/ $1:~/Documents/Projects/focal_adhesions/results/S178A;

echo "Syncing Lin Region Variation Results to $1";
rsync --progress --exclude=**template** -a ../../results/lin_region_variation/ $1:~/Documents/Projects/focal_adhesions/results/lin_region_variation/;

echo "Syncing Simulation Results to $1";
rsync --progress -a ../../results/simulation/ $1:~/Documents/Projects/focal_adhesions/results/simulation/;
