#!/usr/bin/env bash
echo "Syncing Data from Balder";
rsync -a balder:/Volumes/Data/projects/focal_adhesions/data/ ../../data/;
echo "Syncing Results from Balder";
rsync -a balder:/Volumes/Data/projects/focal_adhesions/results/ ../../results/;
