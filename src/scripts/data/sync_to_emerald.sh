#!/usr/bin/env bash
echo "Syncing Directories to Balder";
rsync --exclude-from .balder_exclude -a ../../data emerald:/afs/isis.unc.edu/home/m/b/mbergins/netscr/focal_adhesions/data;
