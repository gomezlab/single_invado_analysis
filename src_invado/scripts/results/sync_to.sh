#!/usr/bin/env bash
echo "Syncing Results Folder to $1";
rsync --progress -a ../../results/* $1:~/Documents/Projects/focal_adhesions/results/;
