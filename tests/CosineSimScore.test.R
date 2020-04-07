# CosineSimScore.test.R
#
# A test of 'CosineSimScore' function in versions above 1.6.0 against the output 
# from the 'CosineSimScore' function in version 1.6.0 as the expected value.

## Load the local source files (code under test)
source('~/Repo/CellScore/R/CosineSimScore.R')
source('~/Repo/CellScore/R/utils.R')

## Load the dependencies
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

   ## Generate cosine similarity score for the combined data with 'CosineSimScore' 
   ## function in the current version and the same function in version 1.6.0 
   ##
   ## NOTE: May take 1-2 minutes on the full eset objects
   ## so we subset it for 4 cell types
   pdata <- pData(eset)
   sel.samples <- pdata$general_cell_type %in% c("ESC", "EC", "FIB", "KER")
   eset.sub <- eset[, sel.samples]

   csExpected <- CellScore::CosineSimScore(eset.sub, cell.change, iqr.cutoff=0.1)
   cs <- CosineSimScore(eset.sub, cell.change, iqr.cutoff=0.1)

    ## Test if the outputs from the 'CosineSimScore' functions in the different 
    ## versions are consistent
   .stopIfNotIdentical(csExpected[[1]], 
      as(cs[[1]], "data.frame"), "csExpected[[1]]", "cs[[1]]", FALSE)
   .stopIfNotIdentical(csExpected[2:5], cs[2:5], "csExpected", "cs", FALSE)
}
