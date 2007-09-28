#!/usr/bin/env bash
echo "Syncing Data to Balder";
rsync -a ../../data/ balder:/Volumes/Data/projects/focal_adhesions/data/;
echo "Syncing Results to Balder";
rsync -a ../../results/ balder:/Volumes/Data/projects/focal_adhesions/results/;
