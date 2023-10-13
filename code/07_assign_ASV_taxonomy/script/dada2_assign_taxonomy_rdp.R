
# assign taxonomy to DADA2 ASVs

#### Libraries ####

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.17")
library(dada2); packageVersion("dada2")

#### direct output ####

outDir = "../output/rdp"
suppressWarnings(dir.create(outDir, recursive = T))

#### find ASV file ####

inputModule = dir("../..", pattern="DADA2_part2", full.names = T)
asvFile = dir(inputModule, pattern="asv-id-sequence.txt", full.names = T, recursive = T)

asv.id.table = read.delim2(asvFile)
message("Read in ", nrow(asv.id.table), " ASVs.")

# assign taxa
taxa <- assignTaxonomy(seqs=asv.id.table$ASV, "~/tax/dada2/rdp_18/rdp_train_set_18.fa.gz", multithread=TRUE)

# Add species.
taxaS <- addSpecies(taxa, "~/tax/dada2/rdp_18/rdp_species_assignment_18.fa.gz")

# if the call to addSpecies() hits error: Error: vector memory exhausted (limit reached?)
# see solutions here: https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached
# or (same solution) here: 
# in short, add line "R_MAX_VSIZE=14Gb " to file ~/.Renviron

taxaS2 = data.frame(taxaS)
taxaS2$ASV = row.names(taxaS)

asvTaxa = merge(asv.id.table, taxaS2, by="ASV")

# match the original order
row.names(asvTaxa) = asvTaxa$asvID
asvTaxa = asvTaxa[asv.id.table$asvID,]

# Save taxa info.
write.table(asvTaxa, 
            file=file.path(outDir, "asvTaxaTable_rdp-v18.txt"), 
            sep="\t", quote=F, row.names = F)

# summary
lines = c(paste0("Started with ", nrow(asv.id.table), " ASVs."),
          paste0("Ended with a taxa table with ", nrow(taxaS2), " rows."))
levelList = c()
for (level in names(taxaS2)){
  levelList = c(levelList, level)
  lines = c(lines, paste0("=============================", level, ":"))
  lines = c(lines, paste0("Number of ASVs with a non-NA value: ", sum(!is.na(taxaS2[,level]))))
  lines = c(lines, paste0("Number of ASVs with a unique taxonomic string: ", length(unique( apply(taxaS2[,levelList, drop=F], 1, paste, collapse=":") ))  ))
  lines = c(lines, paste0("(taxonomic string is ", paste0(levelList, collapse = ":"), ")"))
}

writeLines(lines, file.path(outDir, "summary_rdp-18.txt"))

#### END ####

sessionInfo()

message("Done!")
