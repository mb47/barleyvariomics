import sys
import csv
#import re

## Yun-Yu note: 2020/03/09
## This script is used to validate the quality of low MAF SNPs in the variome project.
## The script reads in a line of vcftools GT[GQ] INFO field output, determine which genotype is the minor allele,
## then separate the GQ score of each individual and label them either as 'major' or 'minor'
## The output file is a long tsv file with the format of $[GQ_score]\t$[minor/major], which is then used to plot the GQ score histogram

## Read file in arg[1]
file=sys.argv[1]


with open(file) as file:
	## Read tsv file using tab as delimiter
	rd = csv.reader(file, delimiter="\t")
	## Read file by row (i.e. per SNP site)
	next(rd)
	for row in rd:
		## Read each element (i.e. per individual) in the row
		## First, determine whether ref or alt is the minor allele
		## Concat the list into a temp string and count occurence of each
		strTemp = ','.join(row)
		refCount = strTemp.count("0/0")
		altCount = strTemp.count("1/1")
		## If alt is the minor allele
		if refCount > altCount:
			for i in row:
				## Parse the genotype field so we can make if judgement and extract the GQ information
				parsed = i.split(':')
				if   parsed[0] == "./.":
					break
				elif parsed[0] == "0/0":
					print("major","\t",parsed[3])
				else:
					print("minor","\t",parsed[3])
		## If ref is the minor allel
		elif refCount < altCount:
			for i in row:
				## Parse the genotype field so we can make if judgement and extract the GQ information
				parsed = i.split(':')
				if   parsed[0] == "./.":
					break
				elif parsed[0] == "0/0":
					print("minor","\t",parsed[3])
				else:
					print("major","\t",parsed[3])
		else:
			break


