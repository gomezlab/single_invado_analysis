#!/usr/bin/env bash
rsync --delete -a balder:/Users/Shared/mbergins/data_store/focal_adhesion_data/data/ ../../data/;
rsync --delete -a balder:/Users/Shared/mbergins/data_store/focal_adhesion_data/results/ ../../data/;
