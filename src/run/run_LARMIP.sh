#PBS -N LARMIP
#PBS -q ns
eval "$(/perm/ms/nl/nkaj/mambaforge/bin/conda shell.bash hook)" 
conda activate ECE

python LARMIP EC-Earth3P-HR highres-future r2i1p2f1

