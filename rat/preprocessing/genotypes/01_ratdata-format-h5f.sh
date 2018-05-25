analysisdir=~/LiMMBoAnalysis/rat/preprocessing/genotypes
datadir=~/data/LiMMBo/rat/rawdata/arrayexpress

python $analysisdir/01_ratdata-format-h5f.py \
    --fname_in HS.hdf5 \
    --fname_out HS_rats \
    --dir $datadir
