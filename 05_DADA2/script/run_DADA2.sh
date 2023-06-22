
# run the dada2 Rscript

mkdir ../log

for BATCH in HiSeq1_4 HiSeq1_8 HiSeq2_3 HiSeq2_4 HiSeq2_7 MiSeq1 MiSeq3 otherRun
do
	echo $BATCH
	Rscript ./dada2_filter_learn_denoise.R $BATCH &> ../log/${BATCH}.log
done

./merge_DADA2_tables.R
