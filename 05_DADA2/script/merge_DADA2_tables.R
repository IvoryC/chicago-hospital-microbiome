
# merge multiple seqtab tables

# Find all of the asv tables in output, 
# merge them, and replace the ASV column names with asv ids.
# Save the ASV seqs to an id file.
# Sort the columns by total abundance, and save the table.

library(dada2)

# # Find the tables
# batchOutput = dir("output", full.names = T, include.dirs = T)
# asvFiles = dir("output", pattern="asv-", full.names = T, recursive = T)
# 
# # read tables
# seqtabs = lapply(asvFiles, read.delim2, row.names=1)
# 
# # merge the tables
# # bigtab1 = do.call(seqtabs, merge)
# bagtab1 = mergeSequenceTables(tables = seqtabs)

asvRSDFiles = dir("../temp", pattern="asv-", full.names = T, recursive = T)
bagtab = mergeSequenceTables(tables = asvRSDFiles)

# # sort the tables
# totals = colSums(bigtab1)
# bigtab = bigtab1[,order(totals, decreasing = T)]

# rename the tables
asvID = paste0("ASV.", 1:ncol(bigtab))
asv.id.table = data.frame(asvID=asvID, ASV=colnames(bigtab))
colnames(bigtab) = asvID

# save the data table
asv.out.File = file.path("../output", paste0("asv-counts.txt"))
write.table(x=bigtab, file=asv.out.File, quote=F, sep="\t", row.names = F)
message("ASV counts table was written to: ", asv.out.File)

# save the table linking asv id to sequence
asv.id.File = file.path("../output", paste0("asv-id-sequence.txt"))
write.table(x=asv.id.table, file=asv.id.File, quote=F, sep="\t", row.names = F)
message("ASV ids and full sequences were written to: ", asv.id.File)

# END

message("Done!")

sessionInfo()
