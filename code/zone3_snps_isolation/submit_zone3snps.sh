#$ -cwd
#$ -j yes
#$ -V
#$ -pe smp 16
#$ -t 1-7
#$ -q ln.q

bcfFile=$1
tempOutDir=$2
sampleFile=815accessions.txt

chrom[1]=chr1H
chrom[2]=chr2H
chrom[3]=chr3H
chrom[4]=chr4H
chrom[5]=chr5H
chrom[6]=chr6H
chrom[7]=chr7H
chroms=${chrom[SGE_TASK_ID]}
region[1]=chr1H:151200000-303200000
region[2]=chr2H:256000000-401600000
region[3]=chr3H:224000000-380000000
region[4]=chr4H:163200000-406500000
region[5]=chr5H:168800000-301600000
region[6]=chr6H:199200000-316000000
region[7]=chr7H:254400000-384000000
regions=${region[SGE_TASK_ID]}

conda activate delAllel

#mkdir $tempOutDir/$chroms
#mkdir $tempOutDir/$chroms.bestTree

## 1. Isolate zone 3 snps
#bcftools view --threads $NSLOTS --regions $regions --samples-file $sampleFile --output-type v --output-file $tempOutDir/$chroms/$chroms.zone3.vcf $bcfFile
#bcftools stats --threads $NSLOTS chr4H.zone3.vcf > $tempOutDir/$chroms/$chroms.zone3.stats

## 1.5 prepare files for PCoA for PAST3
#java utils.snps.vcf.ConvertVCFToFlapjackFormat $tempOutDir/$chroms/$chroms.zone3.vcf $tempOutDir/$chroms/$chroms.zone3.pcoa.dat $tempOutDir/$chroms/$chroms.zone3.pcoa.map

## 2. Convert vcf to phylip format
#python vcf2phylip/vcf2phylip.py -i $tempOutDir/$chroms/$chroms.zone3.vcf

## 3. Model test
#modeltest-ng --datatype nt --input $tempOutDir/$chroms/$chroms.zone3.min4.phy --output $tempOutDir/$chroms/$chroms.modeltest -p $NSLOTS

## 4. ML analysis
# Set model according to AICc criterion
model[1]=GTR+G4
model[2]=TVM+G4
model[3]=TVM+G4
model[4]=GTR+G4
model[5]=TVM
model[6]=TVM
model[7]=TVM+G4
models=${model[SGE_TASK_ID]}

### Run ML only for best tree
#programs/raxml-ng_v0.9.0_linux_x86_64/raxml-ng \
#--msa    $tempOutDir/$chroms/$chroms.zone3.min4.phy \
#--model  $models \
#--prefix $tempOutDir/$chroms.bestTree/zone3 --threads $NSLOTS

# Copy the best tree file and relabell them
#cp $tempOutDir/$chroms.bestTree/zone3.raxml.bestTree $tempOutDir/$chroms.bestTree/$chroms.relablled.bestTree
#project-variome/script/changeBranchLabel.sh $tempOutDir/$chroms.bestTree/$chroms.relablled.bestTree

# Run TreeCluster to see the clustering structure
#TreeCluster.py --input $tempOutDir/$chroms.bestTree/zone3.raxml.bestTree --threshold 0.045 >  $tempOutDir/$chroms.bestTree/$chroms.cluster.txt

# Run bootstraping
qsub submit_bootstrap.sh $chroms $models $tempOutDir

# Plot group by geo-location
#TO DO TO DO TO DO

### Run ML tree with 200 bootstraps
#programs/raxml-ng_v0.9.0_linux_x86_64/raxml-ng --all \
#--msa    $tempOutDir/$chroms/$chroms.zone3.min4.phy \
#--model  $models \
#--prefix $tempOutDir/$chroms/zone3.btsp200 --threads $NSLOTS --bs-trees 200

# Copy the best tree file and relabell them
#cp $tempOutDir/$chroms/zone3.btsp200.raxml.bestTree $tempOutDir/$chroms/$chroms.btsp200.relablled.bestTree
#project-variome/script/changeBranchLabel.sh $tempOutDir/$chroms/$chroms.btsp200.relablled.bestTree

conda deactivate
