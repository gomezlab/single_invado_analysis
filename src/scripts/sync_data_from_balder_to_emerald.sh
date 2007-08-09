#!/usr/bin/env bash
rsync --exclude-from=.emerald_exclude --delete -a balder:/Users/Shared/mbergins/data_store/focal_adhesion_data/ ../../data/;
