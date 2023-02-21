#!/bin/bash

src_data_pre="~/bio/datasets/methylation/20_DMA/02_dmp"
dmr_pre="~/bio/datasets/methylation/20_DMA/03_dmr"

mkdir -p ${dmr_pre}/comb-p


comb-p pipeline \
    -c 4 \
    --seed 0.00005 \
    --dist 5000 \
    --acf-dist 1000 \
    --step 200 \
    -p ${dmr_pre}/comb-p/dex \
    ${src_data_pre}/dmps_svs_for_combp.bed
