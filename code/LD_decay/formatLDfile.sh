#!/bin/bash

file=$1
path=$2

awk -v OFS="\t" 'function abs(v) {return v < 0 ? -v : v} {print $1,$2,$3,$4,$5,$6,$7,abs($2-$5)}' $2/$1 > $2/formatted.$1
