#!/usr/bin/env Rscript

require(BSgenome)
source("~/git/rmaize/R/BSgenomeForge.R")
setwd('~/projects/s3/zhoup-genome/Zmays_B73/21_dbs/bsgenome')
forgeBSgenomeDataPkg("seed", dest='.')


