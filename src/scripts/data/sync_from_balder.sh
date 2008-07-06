#!/usr/bin/env bash
echo "Syncing Directories from Balder";
rsync --exclude-from .balder_exclude -a balder:/Volumes/Data/projects/focal_adhesions/ ../../;
