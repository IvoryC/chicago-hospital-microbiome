#!/bin/zsh

export PATH=$PATH:~/iInstalled/sratoolkit.3.0.5-mac64/bin

while read SRR  ; do
        echo $SRR
        fasterq-dump $SRR -v --split-3 --outdir . 
        gzip ${SRR}.fastq
done < SRR_Acc_List.txt
