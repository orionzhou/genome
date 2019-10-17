require(devtools)
load_all('~/git/rmaize')
require(ggpubr)
dirp = '~/projects/genome'
dird = file.path(dirp, 'data')
dirg = dird

infer_intron <- function(td, zero=T) { # sorted by chrom and start
    #{{{
    chrom = td$chrom[1]; srd = td$srd[1]
    obegs = td$start; oends = td$end
    begs = oends[-nrow(td)]; ends = obegs[-1]
    if(zero == F) begs = begs + 1
    tibble(chrom=chrom, start=begs, end=ends, srd=srd)
    #}}}
}

