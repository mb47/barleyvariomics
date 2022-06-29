#$ -cwd
#$ -j yes
#$ -V

conda activate delAllel
vcfFile=zone3.nonsynonymousSnps.vcf

#### cultivar vs wild ####
## Filtering SNPs by Fst values above 0.8 (cultivars to wild)
fstFile=pairwFst.cultivars.spontaneums.weir.fst
grep -v -w "nan" $fstFile > $fstFile.cleaned
fstFile=$fstFile.cleaned
threshold=0.8
mkdir Cultivar_fixedAllelesList
awk -v threshold=$threshold -v OFS="\t" '$1=="chr1H" && $2 >= 151200000 && $2<= 303200000 && $3 >= threshold {print $1,$2,$2,$3}' $fstFile >  Cultivar_fixedAllelesList/pericentromeric.fixedAlleles.bed
awk -v threshold=$threshold -v OFS="\t" '$1=="chr2H" && $2 >= 256000000 && $2<= 401600000 && $3 >= threshold {print $1,$2,$2,$3}' $fstFile >> Cultivar_fixedAllelesList/pericentromeric.fixedAlleles.bed
awk -v threshold=$threshold -v OFS="\t" '$1=="chr3H" && $2 >= 224000000 && $2<= 380000000 && $3 >= threshold {print $1,$2,$2,$3}' $fstFile >> Cultivar_fixedAllelesList/pericentromeric.fixedAlleles.bed
awk -v threshold=$threshold -v OFS="\t" '$1=="chr4H" && $2 >= 163200000 && $2<= 406400000 && $3 >= threshold {print $1,$2,$2,$3}' $fstFile >> Cultivar_fixedAllelesList/pericentromeric.fixedAlleles.bed
awk -v threshold=$threshold -v OFS="\t" '$1=="chr5H" && $2 >= 168800000 && $2<= 301600000 && $3 >= threshold {print $1,$2,$2,$3}' $fstFile >> Cultivar_fixedAllelesList/pericentromeric.fixedAlleles.bed
awk -v threshold=$threshold -v OFS="\t" '$1=="chr6H" && $2 >= 199200000 && $2<= 316000000 && $3 >= threshold {print $1,$2,$2,$3}' $fstFile >> Cultivar_fixedAllelesList/pericentromeric.fixedAlleles.bed
awk -v threshold=$threshold -v OFS="\t" '$1=="chr7H" && $2 >= 254400000 && $2<= 384000000 && $3 >= threshold {print $1,$2,$2,$3}' $fstFile >> Cultivar_fixedAllelesList/pericentromeric.fixedAlleles.bed

## Intersect between the bed file and LOF annotation
bedtools intersect -a $vcfFile -b Cultivar_fixedAllelesList/pericentromeric.fixedAlleles.bed -wa | sort | uniq > Cultivar_fixedAllelesList/potential.delAlleles.vcf

#==============================================================================

#### landrace vs wild ####
## Filtering SNPs by Fst values above 0.8 (landraces to wild)
fstFile=pairwFst.landraces.spontaneums.weir.fst
grep -v	-w "nan" $fstFile > $fstFile.cleaned
fstFile=$fstFile.cleaned
threshold=0.8
mkdir Landrace_fixedAllelesList
awk -v threshold=$threshold -v OFS="\t" '$1=="chr1H" && $2 >= 151200000 && $2<= 303200000 && $3 >= threshold {print $1,$2,$2,$3}' $fstFile >  Landrace_fixedAllelesList/pericentromeric.fixedAlleles.bed
awk -v threshold=$threshold -v OFS="\t" '$1=="chr2H" && $2 >= 256000000 && $2<= 401600000 && $3 >= threshold {print $1,$2,$2,$3}' $fstFile >> Landrace_fixedAllelesList/pericentromeric.fixedAlleles.bed
awk -v threshold=$threshold -v OFS="\t" '$1=="chr3H" && $2 >= 224000000 && $2<= 380000000 && $3 >= threshold {print $1,$2,$2,$3}' $fstFile >> Landrace_fixedAllelesList/pericentromeric.fixedAlleles.bed
awk -v threshold=$threshold -v OFS="\t" '$1=="chr4H" && $2 >= 163200000 && $2<= 406400000 && $3 >= threshold {print $1,$2,$2,$3}' $fstFile >> Landrace_fixedAllelesList/pericentromeric.fixedAlleles.bed
awk -v threshold=$threshold -v OFS="\t" '$1=="chr5H" && $2 >= 168800000 && $2<= 301600000 && $3 >= threshold {print $1,$2,$2,$3}' $fstFile >> Landrace_fixedAllelesList/pericentromeric.fixedAlleles.bed
awk -v threshold=$threshold -v OFS="\t" '$1=="chr6H" && $2 >= 199200000 && $2<= 316000000 && $3 >= threshold {print $1,$2,$2,$3}' $fstFile >> Landrace_fixedAllelesList/pericentromeric.fixedAlleles.bed
awk -v threshold=$threshold -v OFS="\t" '$1=="chr7H" && $2 >= 254400000 && $2<= 384000000 && $3 >= threshold {print $1,$2,$2,$3}' $fstFile >> Landrace_fixedAllelesList/pericentromeric.fixedAlleles.bed

## Intersect between the bed file and LOF annotation
bedtools intersect -a $vcfFile -b Landrace_fixedAllelesList/pericentromeric.fixedAlleles.bed -wa | sort | uniq > Landrace_fixedAllelesList/potential.delAlleles.vcf
