#!/usr/bin/env bash
echo "Destructively Syncing Data from Balder";
rsync --delete -a balder:/Volumes/Data/projects/focal_adhesions/data/ ../../data/;
echo "Destructively Syncing Results from Balder";
rsync --delete -a balder:/Volumes/Data/projects/focal_adhesions/results/ ../../results/;
