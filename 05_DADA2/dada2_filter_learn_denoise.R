
# Filter seqs in preparation for running DADA2

#### Libraries ####

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.17")
library(dada2); packageVersion("dada2")

#### command line args ####

# This expects an argument telling it which subset of the metadata to work on
# so as to only process this "batch"
# example: HiSeq1_8

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  stop("Expected an argument, such as HiSeq1_4")
}else if (length(args) == 1) {
  batchID = args[1]
}else {
  stop("Too many arguments. Expected one, such as HiSeq1_4")
}
message("Processing reads with run prefix: ", batchID)

#### Find metadata ####

metaFile = "../02_review_metadata_PRJEB14474/output/PRJEB14474_SraRunTable_reduced.txt"
meta = read.delim2(metaFile)
# row.names(data) = data$ID
message("Read metadata from: ", metaFile)

# only keep the files in this "batch"
if (batchID=="ALL") {
  thisBatch = rep(TRUE, nrow(meta))
}else{
  thisBatch = meta$run_prefix..exp.==batchID
}
meta = meta[thisBatch,]
message("Processing ", nrow(meta), " sequence files.")

#### Find sequences ####

gzDir = "/Users/ieclabau/BigPublicData/PRJEB14474"
gzDirs = dir(gzDir, pattern = "^ERR14", full.names = T)
gzFiles = unlist(lapply(gzDirs, dir, pattern="fastq.gz", full.names = T))

# gzFiles = dir(gzDir, pattern="fastq.gz", full.names = T)
ids = gsub(".fastq.gz", "", basename(gzFiles))
names(gzFiles) = ids

gzFiles = gzFiles[names(gzFiles) %in% meta$ID]
ids = names(gzFiles)

if (length(gzFiles)==0) stop("There is no overlap between the availble files and the specified batch.")
#TODO: replace this warning with a stop().
if (sum(thisBatch) > length(gzFiles)) warning("Some metadata rows in this batch did not have a corresponding sequence file.")

#### direct output ####

tmpDir = file.path("temp", batchID)
suppressWarnings(dir.create(tmpDir, recursive = T))

outDir = file.path("output", paste0("batch-", batchID))
suppressWarnings(dir.create(outDir, recursive = T))

#### QC ####

# #skip to save time
# 
# pqp = plotQualityProfile(gzFiles)
# ggsave(plot=pqp, file.path(outDir, "qualityProfile.png"), device = "png")
# pqp

#### filter ####

filtFiles = file.path(tmpDir, basename(gzFiles))
names(filtFiles) = ids

truncLen=120
message("Using a truncation length of: ", truncLen)

filterTable <- filterAndTrim(gzFiles, filtFiles, truncLen=truncLen,
                         maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE)
row.names(filterTable) = ids

filterTableFile = file.path(tmpDir, "filterTable.txt")
write.table(filterTable, filterTableFile, quote=F, sep="\t")
message("Out of ", length(ids), " input files, ", 
        sum(filterTable[,"reads.out"] == 0), 
        " had 0 reads pass the filter. See file ", filterTableFile)

passedFilter = row.names(filterTable)[filterTable[,"reads.out"] > 0]
gzFiles = gzFiles[passedFilter]
filtFiles = filtFiles[passedFilter]
filterTable = filterTable[passedFilter,]
ids = passedFilter
message("Proceding with ", length(passedFilter), " files.")

#### learn errors ####

errorModel <- learnErrors(filtFiles, multithread=TRUE)

# # skip plot to save time
# pe = plotErrors(errA, nominalQ=TRUE)
# ggsave(file.path(outDir,  "plotErrors.png"), plot=pe)
# pe


errorModelFile = file.path(outDir, "dada2ErrorModel.Rdata")
save(errorModel, file=errorModelFile)
message("Learned errors are saved to: ", errorModelFile)

#### assign ASV ####

dadaAs <- dada(filtFiles, err=errorModel, multithread=TRUE)

seqtab = makeSequenceTable(dadaAs)

#### remove chimeras ####

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

perChim = round(1 - sum(seqtab.nochim)/sum(seqtab), 1)
message("About ", perChim * 100, "% of the sequences were the \"chimera\" ASVs.")

message("Moving forward with ", ncol(seqtab.nochim), " ASVs.")

#### filter scarce ASV ####

# an ASV have **more than** 1 read in 1% of samples or an equivalent number of total reads.
# using ceiling and requiring > not just >= means that a count of 1 in 1 in 
# a data set with <100 reads, is still removed.
minCount = ceiling(nrow(seqtab.nochim) * .01)
keepASV = colSums(seqtab.nochim) > minCount
seqtab.filtered = seqtab.nochim[,keepASV]

message("Dropped ", ncol(seqtab.nochim) - ncol(seqtab.filtered), " scarce ASVs.")

message("All ASVs: ", ncol(seqtab))
message("non-chimeric ASVs: ", ncol(seqtab.nochim))
message("ASVs after 1% filter (proceeding): ", ncol(seqtab.filtered))
message("ASVs removed by 1% filter (proceeding): ", sum(colSums(seqtab.nochim) < ceiling(nrow(seqtab.nochim) * .01) ))
message("ASVs removed by 2% filter (hypothetical): ", sum(colSums(seqtab.nochim) < ceiling(nrow(seqtab.nochim) * .02) ))
message("ASVs removed by 3% filter (hypothetical): ", sum(colSums(seqtab.nochim) < ceiling(nrow(seqtab.nochim) * .03) ))
message("ASVs removed by 4% filter (hypothetical): ", sum(colSums(seqtab.nochim) < ceiling(nrow(seqtab.nochim) * .04) ))
message("ASVs removed by 5% filter (hypothetical): ", sum(colSums(seqtab.nochim) < ceiling(nrow(seqtab.nochim) * .05) ))


# png(file.path(outDir, "asv-total-abundance.png"))
# 
# plot(1:ncol(seqtab), log(colSums(seqtab)), las=1, 
#      main=paste0("ASV totals (", batchID, ")"), 
#      ylab="log10 of sum of ASV total counts accross all samples", 
#      xlab="ASVs ordered by abundance")
# mtext(text=paste0(nrow(seqtab), " sequence files"), adj=0)
# mtext(text=paste0(ncol(seqtab), " ASVs"), adj=1)
# abline(h=1:20, col="gray90")
# linesAtX=seq(1000,ncol(seqtab), 1000)
# segments(x0=linesAtX, x1=linesAtX, y0=0, y1=log(colSums(seqtab))[linesAtX], col="gray90")
# points(1:ncol(seqtab), log(colSums(seqtab)))
# 
# dev.off()



#### change ASV to ids ####

if (batchID=="ALL") {
  
  asvID = paste0("ASV.", 1:ncol(seqtab.filtered))
  asv.id.table = data.frame(asvID=asvID, ASV=colnames(seqtab.filtered))
  colnames(seqtab.filtered) = asvID
  
  asv.id.File = file.path(outDir, paste0("asv-id-seq", batchID, ".txt"))
  write.table(x=asv.id.table, file=asv.id.File, quote=F, sep="\t", row.names = F)
  message("ASV ids and full sequences were written to: ", asv.id.File)
  
}

#### save output ####

asvFile = file.path(outDir, paste0("asv-", batchID, ".txt"))
write.table(x=cbind(ID=row.names(seqtab.filtered), seqtab.filtered), 
            file=asvFile, quote=F, sep="\t", row.names = F)



#### tracking ####

track <- cbind(reads.raw = filterTable[,"reads.in"], 
               reads.filter.removed = filterTable[,"reads.in"] - filterTable[,"reads.out"],
               reads.afterFilter = filterTable[,"reads.out"],
               reads.RemovedByDADA2 = filterTable[,"reads.out"] - sapply(dadaAs, function(x) sum(getUniques(x))),
               reads.Counted = sapply(dadaAs, function(x) sum(getUniques(x))), 
               uniqueSeqs=sapply(dadaAs, function(x) length(x$pval)), 
               nASV=apply(seqtab, 1, function(x) sum(x > 0)),
               chimericASVs = sapply(dadaAs, function(dd2) sum( ! dd2$sequence %in% colnames(seqtab.nochim))),
               nASV.nochim = apply(seqtab.nochim, 1, function(x) sum(x > 0)),
               reads.counted.nochim = rowSums(seqtab.nochim),
               scarceASVs = apply(seqtab.filtered, 1, function(x) sum(x > 0)) - apply(seqtab.nochim, 1, function(x) sum(x > 0)),
               nASV.final = apply(seqtab.filtered, 1, function(x) sum(x > 0)),
               reads.counted.final = rowSums(seqtab.filtered))
head(track)
trackFile = file.path(outDir, "trackReadCounts.txt")
write.table(cbind(ID=ids, track), 
            file=trackFile, 
            sep="\t", quote=F, row.names = F)

message("For read count tracking per sample, see file: ", trackFile)

#### END ####

sessionInfo()

message("Done!")
