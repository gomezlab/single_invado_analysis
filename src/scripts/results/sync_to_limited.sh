#!/usr/bin/env bash
echo "Syncing Focal Adhesion to $1";
rsync --progress -a ../../results/focal_adhesions/ $1:~/Documents/Projects/focal_adhesions/results/focal_adhesions/;
echo "Syncing Simulation to $1";
rsync --progress -a ../../results/simulation/ $1:~/Documents/Projects/focal_adhesions/results/simulation/;
