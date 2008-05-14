#!/usr/bin/env bash
echo "Destructively Syncing Directories to Balder";
rsync --exclude-from .balder_exclude --delete -a ../../ balder:/Volumes/Data/projects/focal_adhesions/;
