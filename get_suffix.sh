#!/usr/bin/bash

file=`ls *.gz|head -1`
tmp=${file%.*}
echo ${tmp#*.}