#!/usr/bin/env bash
echo "Syncing Directories from $1";
rsync --exclude-from .balder_exclude -a $1:~/Documents/Projects/focal_adhesions/ ../../;
