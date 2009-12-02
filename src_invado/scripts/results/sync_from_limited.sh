#!/usr/bin/env bash
echo "Syncing Focal Adhesion from $1";
rsync --progress -a $1:~/Documents/Projects/focal_adhesions/results/* ../../results/;
