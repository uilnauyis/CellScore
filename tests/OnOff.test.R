# OnOff.test.R
#
# A test of 'OnOff' function in versions above version 1.6.0 against the output 
# from the 'OnOff' function in the version 1.6.0 as the expected value

## Load the source files (code under test)
source('~/Repo/CellScore/R/OnOff.R')
source('~/Repo/CellScore/R/utils.R')

## Load the dependences
library(Biobase)
library(hgu133plus2CellScore) # eset.std
library(SummarizedExperiment)

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

    ## Generate two marker lists with 'OnOff' function in the current version and 
    ## the same function in version 1.6.0 
    group.OnOff <- OnOff(eset, cell.change, out.put="marker.list")
    group.OnOffExpected <- CellScore::OnOff(eset, cell.change, out.put="marker.list")

    ## Calculate the on/off scores for individual samples with with the 'OnOff' 
    ## function in the current version and the same function in version 1.6.0    
    individ.OnOff <- OnOff(eset, cell.change, out.put="individual")
    individ.OnOffExpected <- CellScore::OnOff(eset, cell.change, out.put="individual")

    ## Test if the outputs from the 'OnOff' functions in the different versions are
    ## consistent
    .stopIfNotIdentical(group.OnOff, group.OnOffExpected, "group.OnOff",
                        "group.OnOffExpected", FALSE)
    .stopIfNotIdentical(individ.OnOff, individ.OnOffExpected, "iindivid.OnOff",
                       "individ.OnOffExpected", FALSE)
}


