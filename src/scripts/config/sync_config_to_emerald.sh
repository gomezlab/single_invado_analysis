#!/usr/bin/env bash
echo "Syncing Config to Balder";
rsync -a ../../data/config/ emerald:~/netscr/focal_adhesions/data/config/;
