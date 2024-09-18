##SET UP -> required packages


if (!require("foreach", quietly = TRUE))
  install.packages("foreach");

if (!require("XML", quietly = TRUE))
  install.packages("XML")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages('devtools')

BiocManager::install("Rsamtools")
BiocManager::install("rtracklayer")

BiocManager::install("BSgenome")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

BiocManager::install("affxparser")


install.packages("remotes")
remotes::install_github("gustaveroussy/apt.oncoscan.2.4.0")

#https://rdrr.io/github/gustaveroussy/EaCoN/f/README.md

install.packages("https://zenodo.org/record/5494853/files/OncoScanCNV.na33.r2_0.1.0.tar.gz", repos = NULL, type = "source")


# for some reason dl always breaks
#wget https://zenodo.org/record/5494853/files/OncoScanCNV.na33.r2_0.1.0.tar.gz
#R -e "install.packages('./OncoScanCNV.na33.r2_0.1.0.tar.gz', repos = NULL, type = 'source')"


BiocManager::install("rhdf5")
BiocManager::install("dplyr")

install.packages("RColorBrewer")
install.packages("tidyr")