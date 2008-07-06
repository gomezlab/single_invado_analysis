#!/usr/bin/env bash
echo "Syncing Config from Balder";
rsync -a balder:/Volumes/Data/projects/focal_adhesions/data/config ../../data/config;
