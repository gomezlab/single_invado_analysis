#!/usr/bin/env bash
echo "Syncing Config from $1";
rsync -a $1:~/Documents/Projects/focal_adhesions/data/config/ ../../data/config;
