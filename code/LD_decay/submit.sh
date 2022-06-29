#$ -cwd
#$ -j yes
#$ -V

./formatLDfile.sh r2.cultivars.ld cultivars
./formatLDfile.sh r2.landraces.ld landraces
./formatLDfile.sh r2.spontaneums.ld spontaneums
