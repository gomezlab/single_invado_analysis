#!/usr/bin/env bash
echo "Syncing Data to Balder";
rsync --exclude-from .balder_exclude -a ../../data/ balder:/Volumes/Data/projects/focal_adhesions/data/;
echo "Syncing Results to Balder";
rsync --exclude-from .balder_exclude -a ../../results/ balder:/Volumes/Data/projects/focal_adhesions/results/;
