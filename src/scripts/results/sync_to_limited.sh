#!/usr/bin/env bash
echo "Syncing Focal Adhesion Data to $1";
rsync --progress --exclude=*visualizations* -a ../../results/focal_adhesions/ $1:~/Documents/Projects/focal_adhesions/results/focal_adhesions/;

echo "Syncing Focal Adhesion Visualizations to $1";
rsync --progress --include=**visualizations** --include=*/ --exclude=** --prune-empty-dirs -a ../../results/focal_adhesions/ $1:~/Documents/Projects/focal_adhesions/results/focal_adhesions/;

echo "Syncing Simulation to $1";
rsync --progress -a ../../results/simulation/ $1:~/Documents/Projects/focal_adhesions/results/simulation/;
