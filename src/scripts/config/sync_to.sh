#!/usr/bin/env bash
echo "Syncing Config to $1";
rsync -a ../../data/config/ $1:~/Documents/Projects/focal_adhesions/data/config;
