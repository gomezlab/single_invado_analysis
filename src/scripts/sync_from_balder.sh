#!/usr/bin/env bash
echo "Syncing Data from Balder";
rsync --exclude-from .balder_exclude -a balder:/Volumes/Data/projects/focal_adhesions/data/ ../../data/;
echo "Syncing Results from Balder";
rsync --exclude-from .balder_exclude -a balder:/Volumes/Data/projects/focal_adhesions/results/ ../../results/;
