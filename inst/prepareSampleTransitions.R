#' Prepare sample transitions
#'
#' @export
#' @examples
#' sampleTransitions <- prepareSampleTransitions()
prepareSampleTransitions <- function() {

    dat <- read.csv(system.file("extdata",
                                  "sampleMetadata.csv",
                                  package = "CellScore"),
                    header = TRUE)

    uniqueCellTypes <- unique(c(as.array(dat$general_cell_type),
                                as.array(dat$target_cell_type),
                                as.array(dat$parental_cell_type)))

    print(uniqueCellTypes)

    cellTypeMap <- getCellTypeMap

    uniqueTransitions <- unique(dat[dat$category == 'test',
                                    c('parental_cell_type',
                                      'general_cell_type',
                                      'target_cell_type')])

    colnames(uniqueTransitions) <- c("start",	"test", "target")

    ## Exclude samples without parental cell types or without target cell type
    uniqueTransitions <- uniqueTransitions[!is.na(uniqueTransitions$start) &
                                             !is.na(uniqueTransitions$target) &
                                             uniqueTransitions$start != '' &
                                             uniqueTransitions$target != '', ]

    ## Replace full cell type names with abbreviations
    uniqueTransitionsAbbrv <- as.data.frame(lapply(uniqueTransitions,
                                                   function(element) {
                                                      (cellTypeMap[element])
                                                   }))


    tempDirPath <- paste(getwd(), 'temp', sep = '/')
    dir.create(tempDirPath)

    write.table(uniqueTransitionsAbbrv,
                file = paste(getwd(),
                             'temp',
                             'exampleTransitionsHsapienData.tsv',
                             sep = '/'),
                quote = FALSE,
                sep = '\t')

    uniqueTransitionsAbbrv
}

getCellTypeMap <- function() {
    ## Hardcode abbreviations for each cell type
    cellTypes <- c(
        "RPE",
        "ESC",
        "ESC.d",
        "ESC.t",
        "iPS",  ## check the comment
        "IPS",
        "URI",
        "PBMC",
        "HUVEC",

        "derived cardiomyocyte",
        "hepatocyte-like cell",
        "hepatocyte",
        "fetal liver",
        "derived dopaminergic neuron",
        "derived neuron",
        "spinal cord",
        "derived neural progenitor cell",
        "derived RPE",
        "fetal RPE",
        "fetal heart",
        "derived neural stem cell",
        "naive IPS",
        "naive ESC",
        "fetal brain",
        "derived pancreatic b-like cell",
        "derived endoderm",
        "derived foregut endoderm",
        "derived hepatic endoderm",
        "derived early hepatocyte",
        "derived late hepatocyte",
        "fibroblast",
        "induced neuron",
        "fibroblast.treated",
        "cardiomyocytes",
        "hepatocytes",
        "dopaminergic neurons",
        "neurons",
        "neual progenitor cell",
        "neural stem cells",
        "pancreatic beta cell",
        "endoderm",
        "foregut endoderm",
        "hepatic endoderm",
        "neuron",
        "lymphocyte")

    abbrv <- c(
        ## Original abbreviations
        "RPE",
        "ESC",
        "ESC.d",
        "ESC.t",
        "IPS",
        "IPS",
        "URI",
        "PBMC",
        "HUVEC",
        ## Abbreviations named by me
        "dCDM",
        "HPTL",
        "HPT",
        "fLIV",
        "dDPNEU",
        "dNEU",
        "SPC",
        "dNEUPR",
        "dRPE",
        "fRPE",
        "fHRT",
        "dNEUST",
        "nIPS",
        "nESC",
        "fBRN",
        "dPCB",
        "dEDD",
        "dFEDD",
        "dHPTEDD",
        "dEHPT",
        "dLHPT",
        "FIB",
        "INEU",
        "FIB.t",
        "CDM",
        "HPT",
        "DPNEU",
        "NEU",
        "NEUPR",
        "NEUST",
        "PCB",
        "EDD",
        "FEDD",
        "HPTEDD",
        "NEU",
        "LPC")

    cellTypeMap <- setNames(abbrv, cellTypes)
    cellTypeMap
}
