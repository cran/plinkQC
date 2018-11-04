## ----setup knitr, include = FALSE----------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----check ancestry, eval=TRUE, fig.height=3, fig.width=5, fig.align='center'----
library(plinkQC)
package.dir <- find.package('plinkQC')
indir <- file.path(package.dir,'extdata')
name <- 'data'
refname <- 'HapMapIII'
prefixMergedDataset <- paste(name, ".", refname, sep="")

exclude_ancestry <-
    evaluate_check_ancestry(indir=indir, name=name,
                            prefixMergedDataset=prefixMergedDataset,
                            refSamplesFile=paste(indir, "/HapMap_ID2Pop.txt",
                                                 sep=""), 
                            refColorsFile=paste(indir, "/HapMap_PopColors.txt",
                                                sep=""),
                            interactive=TRUE)

