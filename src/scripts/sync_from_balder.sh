#!/usr/bin/env bash
rsync --delete -a balder:/Volumes/Data/projects/focal_adhesions/data/ ../../data/;
rsync --delete -a balder:/Volumes/Data/projects/focal_adhesions/results/ ../../results/;
