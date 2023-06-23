
# Initial look via Phyloseq

library(phyloseq); packageVersion("phyloseq")
## [1] ‘1.44.0’
library(Biostrings); packageVersion("Biostrings")
## [1] ‘2.68.1’
library(ggplot2); packageVersion("ggplot2")
## [1] ‘3.4.2’

#### direct output ####

tmpDir = "../temp"
suppressWarnings(dir.create(tmpDir, recursive = T))

outDir = "../output"
suppressWarnings(dir.create(outDir, recursive = T))

#### metadata ####

metaFile = "../../02_review_metadata_PRJEB14474/output/PRJEB14474_SraRunTable_reduced.txt"
meta = read.delim2(metaFile)
row.names(meta) = meta$ID
# row.names(data) = data$ID
message("Read metadata from: ", metaFile)

#### taxa data ####

taxaFile = "../../06_assign_ASV_taxonomy/output/asvTaxaTable.txt"
taxa = read.delim2(taxaFile, row.names = "asvID")
# drop the ASV column
taxaLim = taxa[,grep("ASV", names(taxa), invert = T, value = T)]
taxaMat = as.matrix(taxaLim)

#### counts ####

countsFile = "../../05_DADA2/output/asv-counts.txt"
counts = read.delim2(countsFile, row.names = 1)
dim(counts)

#### trim down data ####

sampleTypesOfInt = c("blank control", "Hand", "Axilla", "Cold Tap Water", "Hot Tap Water", "Countertop")
meta = meta[meta$SAMPLE_TYPE %in% sampleTypesOfInt,]
counts = counts[meta$ID,]

splitT = split(counts, f=meta$SAMPLE_TYPE)
splitTsort = lapply(splitT, function(df) df[order(rowSums(df), decreasing = T),])
splitTsortLim = lapply(splitTsort, function(df) {n = min(c(nrow(df), 10)); df[1:n,] } )
# put it back together
counts = do.call(rbind, splitTsortLim)
# fix row names
row.names(counts) = sapply(strsplit(row.names(counts), split=".", fixed=T), "[[", 2)

# This is just a sanity check, only keep the 100 most abundance taxa,
counts = counts[,order(colSums(counts), decreasing = T)]
counts = counts[, 1:100]

# only keep rows/columns that have at least one read.
counts = counts[rowSums(counts) > 0,]

dim(counts)

# match up other tables
meta = meta[row.names(counts),]
taxaMat = taxaMat[names(counts),]

#### Phyloseq ####

ps <- phyloseq(otu_table(counts, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxaMat))

# works, but not very valuable.
plot_richness(ps, x="study_day", measures=c("Shannon", "Simpson"), color="SAMPLE_TYPE") + facet_wrap("surcat")

ggsave(file.path(outDir, "richness.png"))

### ordination
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
po = plot_ordination(ps.prop, ord.nmds.bray, color="SAMPLE_TYPE", title="Bray NMDS")
po

ggsave(po, filename = file.path(outDir, "ordination.png"))
# po + geom_text(aes(label=ID), nudge_y = .02, nudge_x = -5)
# wow, sample ERR1459412 really stands out.  ...it only has a few reads.


# barplot
top20 <- names(counts)[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
pb = plot_bar(ps.top20, x="ID", fill="Family") + facet_wrap(~SAMPLE_TYPE, scales="free_x")
pb

ggsave(file.path(outDir, "barplot_family.png"), pb)








#### END ####

sessionInfo()

message("Done!")


