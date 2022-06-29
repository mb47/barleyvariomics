#!/bin/bash

group[1]=cultivars
group[2]=landraces
group[3]=spontaneums

echo -e CHR'\t'LOC'\t'START'\t'END'\t'VAR'\t'SFS'\t'LD'\t'MU'\t'GROUP > combined_report.txt

for i in {1..3};do
	for j in {1..7};do
	chr=chr${j}H
	group=${group[i]}
	tail -n +2 RAiSD_Report.$group.$chr | awk -v chr=$chr -v OFS="\t" -v group=$group '{print chr,$0,group}' >> combined_report.txt
	done
done
