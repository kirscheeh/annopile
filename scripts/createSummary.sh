#!/bin/bash

CONFIG=path.conf

source $CONFIG

while read -ra SAMPLE; do
    IFS="," read -ra SAMPLE_COMP <<< $SAMPLE
    
    BAM=${SAMPLE_COMP[0]}
    SAMPLE_ID=${SAMPLE_COMP[1]}
    
    python3 $REPO_PATH/scripts/pileups.py $RESULTS $BAM $SAMPLE_ID 
    
    docker run --rm -v $PARENT_ANNOVAR:/ref -v $RESULTS:/run -v $REPO_PATH/scripts/:/scripts -e SAMPLE_ID=$SAMPLE_ID pileup_annovar:1.0 bash /scripts/annovar.sh $SAMPLE_ID

    python3 $REPO_PATH/scripts/summary.py $RESULTS/${SAMPLE_ID}.xlsx $RESULTS/${SAMPLE_ID}.hg19_multianno.csv

done < input/SOI.csv