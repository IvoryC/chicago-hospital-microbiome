
# merge multiple seqtab tables

# Find all of the asv tables in output, 
# merge them, and replace the ASV column names with asv ids.
# Save the ASV seqs to an id file.
# Sort the columns by total abundance, and save the table.

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.17")
library(dada2); packageVersion("dada2")

#### direct output ####

outDir = "../output"
suppressWarnings(dir.create(outDir))

#### Find the tables #### 

# get the output of "05_DADA2_part1" but remember the number could change.
inputModule = dir("../..", pattern="DADA2_part1", full.names = T, include.dirs = T)

pathEndings = c(
  "HiSeq1_4/asv-HiSeq1_4.RDS",
  "HiSeq1_8/asv-HiSeq1_8.RDS",
  "HiSeq2_3/asv-HiSeq2_3.RDS",
  "HiSeq2_4/asv-HiSeq2_4.RDS",
  "HiSeq2_7/asv-HiSeq2_7.RDS",
  "MiSeq1/asv-MiSeq1.RDS",
  "MiSeq3/asv-MiSeq3.RDS",
  "otherRun/asv-otherRun.RDS")
asvRSDFiles = file.path(inputModule, "temp", pathEndings)

# read, merge sort tables
bigtabAll = mergeSequenceTables(tables = asvRSDFiles)

in.nASVpersample = apply(bigtabAll, 1, function(x) sum(x > 0))

#### remove chimeras ####

bigtab <- removeBimeraDenovo(bigtabAll, method="consensus", multithread=TRUE, verbose=TRUE)

percentChim = round(1 - sum(bigtab)/sum(bigtabAll), 1)
message("About ", percentChim * 100, "% of the sequences were the \"chimera\" ASVs.")

message("Moving forward with ", ncol(bigtab), " ASVs.")

nASVpersample = apply(bigtab, 1, function(x) sum(x > 0))

####  rename the tables #### 
asvID = paste0("ASV.", 1:ncol(bigtab))
asv.id.table = data.frame(asvID=asvID, ASV=colnames(bigtab))
colnames(bigtab) = asvID

####  save the data table #### 
asv.out.File = file.path(outDir, paste0("asv-counts.txt"))
write.table(x=cbind(ID=row.names(bigtab), bigtab), 
            file=asv.out.File, quote=F, sep="\t", row.names = F)
message("ASV counts table with ", nrow(bigtab), " samples and ", ncol(bigtab), " ASV ids was written to: ", asv.out.File)

saveRDS(bigtab, file.path(outDir, paste0("asv-counts.RDS")))

####  save the table linking asv id to sequence #### 
asv.id.File = file.path(outDir, paste0("asv-id-sequence.txt"))
write.table(x=asv.id.table, file=asv.id.File, quote=F, sep="\t", row.names = F)
message("ASV ids and full sequences were written to: ", asv.id.File)


#### tracking #### 

# all of these tables have an identical format for the columns, and independent rows.
trackPathEndings = c(
  "batch-HiSeq1_4/trackReadCounts.txt",
  "batch-HiSeq1_8/trackReadCounts.txt",
  "batch-HiSeq2_3/trackReadCounts.txt",
  "batch-HiSeq2_4/trackReadCounts.txt",
  "batch-HiSeq2_7/trackReadCounts.txt",
  "batch-MiSeq1/trackReadCounts.txt",
  "batch-MiSeq3/trackReadCounts.txt",
  "batch-otherRun/trackReadCounts.txt")
trackFiles = file.path(inputModule, "output", trackPathEndings)

trackMulti = lapply(trackFiles, read.delim2)
track = do.call(rbind, trackMulti)

track$chimericASVs = nASVpersample - in.nASVpersample
track$nASV.part2 = nASVpersample
track$reads.part2 = rowSums(bigtab)

# save the tracking data
trackFile = file.path(outDir, paste0("trackReadCounts.txt"))
write.table(x=track, file=trackFile, quote=F, sep="\t", row.names = F)
message("Updated read count tracking table was saved to: ", trackFile)

# END

message("Done!")

sessionInfo()
