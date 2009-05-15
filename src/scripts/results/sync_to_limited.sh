#!/usr/bin/env bash
echo "Syncing Focal Adhesion Data to $1";
rsync --progress --exclude=**visualizations** -a ../../results/focal_adhesions/ $1:~/Documents/Projects/focal_adhesions/results/focal_adhesions/;

rsync --progress -a ../../results/S178A/ $1:~/Documents/Projects/focal_adhesions/results/S178A;

rsync --progress --exclude=**template** -a ../../results/lin_region_variation/ $1:~/Documents/Projects/focal_adhesions/results/lin_region_variation/;

echo "Syncing Focal Adhesion Visualizations to $1";
rsync --progress --include=**visualizations** --include=*/ --exclude=** --prune-empty-dirs -a ../../results/focal_adhesions/ $1:~/Documents/Projects/focal_adhesions/results/focal_adhesions/;

echo "Syncing Simulation to $1";
rsync --progress -a ../../results/simulation/ $1:~/Documents/Projects/focal_adhesions/results/simulation/;
