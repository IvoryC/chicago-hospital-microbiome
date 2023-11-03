# Self find (comparative self-discrimination) (score) functions
# In a meeting on Wednesday Nov 1, Anthony Fodor, Hyotae, Li Ma and Ivory (myself)
# discussed weather or not Decontam was making any difference in the 
# Chicago Hospital dataset.
# 

#### libraries ####
tellme <- function(name){print(paste0("Package ", name, " version: ", packageVersion(name)))}

library(tidyr); tellme("tidyr")
suppressPackageStartupMessages(library(dplyr)); tellme("dplyr")
library(ggplot2); tellme("ggplot2")

#### main functions ####


selfSelectionScore <- function(pairs, counts, title="", showPlots = T, method = "pearson"){
  # pairs - a dataframe or matrix with two columns and two or more rows.
  #.        The two sample ids in each row should have some relationship (such as
  #.        being from the same patient).  There should be no duplicates across rows.
  #.        Some functions assume that row.names are set to the first column.
  # counts - counts per sample, with features as columns and samples as rows.
  #.        The row names of counts should include all sample ids in both columns
  #.        of pairs.
  # showPlots - should plots be printed when the function is run
  # method - passed to cor() via corrMatrix
  corMat = corrMatrix(pairs, counts, method=method)
  
  res = calcRes(corMat, pairs)
  
  heatmap = correlationTableToHeatMap(corMat, title=title)
  points = correlationPoints(corMat, title=title)
  
  if (showPlots){
    show(heatmap)
    show(points)
  }
  
  return(list(correlatoins=corMat, table=res, heatmap=heatmap, points=points))
}



corrMatrix <- function(pairs, counts, method = "pearson"){
  # pairs - a dataframe or matrix with two columns and two or more rows.
  #.        The two sample ids in each row should have some relationship (such as
  #.        being from the same patient).  There should be no duplicates across rows.
  #.        Some functions assume that row.names are set to the first column.
  # counts - counts per sample, with features as columns and samples as rows.
  #.        The row names of counts should include all sample ids in both columns
  #.        of pairs.
  # [value] a matrix with rows corresponding to the first pairs column (left), and columns
  #.        corresponding to the second column (right). Values in the table are correlations.
  left = pairs[,1]
  right = pairs[,2]
  allSamples = unlist(pairs)
  corMatrix = cor(t(counts[allSamples,]), use="pairwise.complete.obs", method = method)
  # cor.co = as.data.frame(cor.co)
  corMat = corMatrix[left, right]
  return(corMat)
}

calcRes <- function(corMat, pairs){
  if (is.null(pairs)){
    pairs = extractPairs(corMat)
  }
  res = data.frame(left=NA, matchCor=0L, unmatchMedian=0L)
  for (left in row.names(pairs)){
    right = pairs[left, 2]
    matchCor = corMat[left, right]
    others = setdiff(pairs[,2], right)
    unmatchMedian = median(unlist(corMat[left, others]))
    res = rbind(res, c(left, matchCor, unmatchMedian))
  }
  # remove starter row
  res = filter(res, !is.na(left))
  # calculate scores
  res$unmatchMedian = as.numeric(res$unmatchMedian)
  res$matchCor = as.numeric(res$matchCor)
  res = mutate(res, score = matchCor / unmatchMedian)
}

extractPairs <- function(corMat){
  # table - a matrix with rownames and colnames in order such that pairs 
  #.        are in the same order on both sides.
  # table = as.data.frame(table)
  pairs = data.frame(left=row.names(corMat),
                     right=colnames(corMat))
  row.names(pairs) = pairs[,1]
  return(pairs)
}


#### plot functions ####

correlationTableToHeatMap <- function(cor.table,
                                     pdfFilename=NULL, 
                                     title="correlations", 
                                     midpointVal=.5, blueVal=0, redVal=1, 
                                     xlab="right", ylab="left"){
  # cor.table - a matrix with row names and column names, and numeric values.
  # title - title for the plot
  # pdfFilename - if supplied, image is saved to a file.
  # blueVal - values worse than this will all be shown in blue.
  # redVal - values higher than this will be shown in red.
  # midpointVal - these will be shown in white, as a middle point between the blue and red gradients.
  # all p-values better than this will be shown in solid red.

  # max out at limits
  cor.table[cor.table > redVal] = redVal
  cor.table[cor.table < blueVal] = blueVal
  longTab = cor.table %>% 
    as.data.frame() %>%
    mutate(left = row.names(cor.table) ) %>%
    pivot_longer(cols=-left, names_to = "right", values_to="correlation")
  longTab$right = factor(x=as.character(longTab$right), levels=colnames(cor.table))
  longTab$left = factor(x=as.character(longTab$left), levels=row.names(cor.table))
  gp = ggplot(longTab, aes(x=right, y=left, fill=correlation)) +
    geom_tile(col="gray") +
    scale_fill_gradient2(low = "white", high = "red", mid = "yellow",
                         midpoint = midpointVal,
                         limit = c(blueVal,redVal),
                         space = "Lab",
                         name="correlation",
                         na.value = gray(.9)) +
    theme_minimal()+ # minimal theme
    ggtitle(title) +
    xlab(xlab) +
    ylab(ylab) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 8, hjust = 1)) +
    coord_fixed()
  
  if (!is.null(pdfFilename)){
    # pdf(pdfFilename, width=30, height = 25)
    pdf(pdfFilename, 
        height = 3 + nrow(cor.table) / 8,
        width = 5 + ncol(cor.table) / 8)
    show(gp)
    dev.off()
    message("See file: ", pdfFilename)
  }
  
  return(gp)
}





correlationPoints <- function(corMat, pairs=NULL, res=NULL, title=""){
  if (is.null(pairs)){
    pairs = extractPairs(corMat)
  }
  if (is.null(res)){
    res=calcRes(corMat, pairs)
  }
  
  longTab = corMat %>% 
    as.data.frame() %>%
    mutate(left= row.names(corMat) ) %>%
    pivot_longer(cols=-left, names_to = "right", values_to="correlation") %>%
    mutate(isMatch = ifelse(pairs[left,2]==right,"yes", "no"))
  
  points = ggplot(longTab, aes(x=left, y=correlation, col=isMatch)) +
    scale_color_manual(values=c(yes="red", no="gray")) +
    geom_jitter(width = .1, height = 0) + 
    theme_minimal() +
    geom_hline(yintercept = c(0,1), col=gray(.5)) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    ylim(0,1) +
    ggtitle(title) + 
    annotate("point", x=res$left, y=res$unmatchMedian, shape="-", size=8) +
    annotate("text", x=res$left, label=signif(res$score,2), 
             y=0, vjust=1.1, hjust=0)
  return(points)
}



#### other handy functions ####

lognorm <- function(table, log10.do = TRUE){
  # table - a table with rows for samples and columns for features
  #         samples names are row names.
  #         feature names are column names.
  #         the table values are all numeric
  sampleSums = rowSums(table)
  meanSampleDepth = mean(sampleSums)
  sampleProportion = t( apply(table, 1, function(row) row / sum(row)) )
  t2 = sampleProportion * meanSampleDepth
  if (log10.do) t3 = log10( t2 + 1 )
  else t3 = t2 + 1
  return( t3 )
}