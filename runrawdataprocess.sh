#!/bin/bash

samplelist="../filereport_read_run_PRJNA348749_tsv.txt"

cat $samplelist | cut -f 4 | while read samp
do
bash rawdataprocess.parseagr.sh -s $samp -1 ../rawdata/${samp}_1.fastq.gz -2 ../rawdata/${samp}_2.fastq.gz -o ./ -f true -t 4

done
