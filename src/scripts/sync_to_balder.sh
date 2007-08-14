#!/usr/bin/env bash
rsync --delete -a ../../data/ balder:/Users/Shared/mbergins/data_store/focal_adhesion_data/data;
rsync --delete -a ../../results/ balder:/Users/Shared/mbergins/data_store/focal_adhesion_data/results;
