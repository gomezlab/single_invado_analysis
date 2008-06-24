#!/usr/bin/env bash
echo "Syncing Config to Balder";
rsync -a ../../data/config balder:/Volumes/Data/projects/focal_adhesions/data/config;
