#!/usr/bin/env bash
echo "Destructively Syncing Directories from $1";
rsync --exclude-from .balder_exclude --delete -a $1:/Volumes/Data/projects/focal_adhesions/ ../../;
