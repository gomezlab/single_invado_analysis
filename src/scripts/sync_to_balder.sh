#!/usr/bin/env bash
rsync --delete -a ../../data/ balder:/Volumes/Data/projects/focal_adhesions/data/;
rsync --delete -a ../../results/ balder:/Volumes/Data/projects/focal_adhesions/results/;
