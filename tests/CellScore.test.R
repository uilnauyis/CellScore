# CellScore.test.R
#
# An integration test of 'CosineSimScore', 'OnOff' and 'CellScore' functions in 
# versions above 1.6.0 against the output from the same functions in version 1.6.0 
# as the expected value.

## Load the local source files (code under test)
source('~/Repo/CellScore/R/OnOff.R')
source('~/Repo/CellScore/R/CosineSimScore.R')
source("~/Repo/CellScore/R/CellScore.R")
source('~/Repo/CellScore/R/utils.R')

## Load the dependences
library(Biobase) 
library(hgu133plus2CellScore) # eset.std
library(SummarizedExperiment)
library(lsa) ## cosine

## Locate the external data files in the CellScore package
rdata.path <- system.file("extdata", "eset48.RData", package = "CellScore")
tsvdata.path <- system.file("extdata", "cell_change_test.tsv",
                             package = "CellScore")

if (file.exists(rdata.path) && file.exists(tsvdata.path)) {

   ## Load the expression set with normalized expressions of 48 test samples
   load(rdata.path)

   ## Import the cell change info for the loaded test samples
   cell.change <- read.delim(file= tsvdata.path, sep="\t",
                             header=TRUE, stringsAsFactors=FALSE)

   ## Combine the standards and the test data
   eset <- combine(eset.std, eset48)

   ## Generate cosine similarity for the combined data with function 
   ## 'CosineSimScore' in current version
   ## NOTE: May take 1-2 minutes on the full eset object
   ## so we subset it for 4 cell types
   pdata <- pData(eset)
   sel.samples <- pdata$general_cell_type %in% c("ESC", "EC", "FIB", "KER")
   eset.sub <- eset[, sel.samples]

   cs <- CosineSimScore(eset.sub, cell.change, iqr.cutoff=0.1)

   ## Generate the on/off scores for the combined data with function 'OnOff' in 
   ## current version
   individ.OnOff <- OnOff(eset.sub, cell.change, out.put="individual")

   ## Generate the CellScore values for all samples with function 'CellScore' in
   ## current version
   cellscore <- CellScore(eset.sub, cell.change, individ.OnOff$scores,
                          cs$cosine.samples)


   ## Generate cosine similarity for the combined data with function 
   ## 'CosineSimScore' in version 1.6.0
   csExpected <- CellScore::CosineSimScore(eset.sub, cell.change, iqr.cutoff=0.1)

   ## Generate the on/off scores for the combined data with function 'OnOff' in 
   ## version 1.6.0
   individ.OnOffExpected <- CellScore::OnOff(eset.sub, cell.change, 
                                             out.put="individual")

   ## Generate the CellScore values for all samples with function 'CellScore' in
   ## version 1.6.0
   cellscoreExpected <- CellScore::CellScore(eset.sub, cell.change, 
                                  individ.OnOffExpected$scores, 
                                  csExpected$cosine.samples)

   ## Test if output from the different versions are consistent
   .stopIfNotIdentical(cellscoreExpected, cellscore, "cellscoreExpected",
                       "cellscore", FALSE)
}
