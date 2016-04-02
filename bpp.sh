#!/bin/bash
#
date
# build-number++
file="build-number"
bn=`cat $file`
echo $[bn+1] > $file

