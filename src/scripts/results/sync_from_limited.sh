#!/usr/bin/env bash
echo "Syncing Focal Adhesion from $1";
rsync -a $1:~/Documents/Projects/focal_adhesions/results/focal_adhesions/ ../../results/focal_adhesions;
echo "Syncing FAK from $1";
rsync -a $1:~/Documents/Projects/focal_adhesions/results/FAK/ ../../results/FAK;
