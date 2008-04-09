#!/usr/bin/env bash
echo "Destructively Syncing Data to Balder";
rsync --exclude-from ../.balder_exclude --delete -a ../../data/ balder:/Volumes/Data/projects/focal_adhesions/data/;
echo "Destructively Syncing Results to Balder";
rsync --exclude-from ../.balder_exclude --delete -a ../../results/ balder:/Volumes/Data/projects/focal_adhesions/results/;
