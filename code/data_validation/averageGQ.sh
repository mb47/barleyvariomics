#$ -cwd
#$ -j yes
#$ -V

conda activate R_env
Rscript plotting.r
