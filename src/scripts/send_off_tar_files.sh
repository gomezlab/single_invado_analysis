#!/usr/bin/env bash

tar zcf - ../../results/focal_adhesions/ | ssh $1 "cat - > fa.tar.gz"
tar zcf - ../../results/S178A/ | ssh $1 "cat - > S178A.tar.gz"
