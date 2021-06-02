require(devtools)
load_all('~/git/rmaize')
require(ggpubr)
dirp = '~/projects/genome'
dird = file.path(dirp, 'data')
dirg = dird
gts3 = c("B73",'Mo17','W22')
gts6 = c("B73",'Mo17','W22','B73xMo17','W22xB73','W22xMo17')
gts25 = str_split("B73 B97 CML322 CML333 CML52 CML69 DK105 EP1 F7 Il14H Ki11 Ki3 M162W M37W Mo17 Mo18W MS71 NC350 NC358 Oh43 Oh7B P39 PH207 Tx303 W22", " ")[[1]]
gts31 = str_split("Mo17 W22 B97 CML52 CML69 CML103 CML228 CML247 CML277 CML322 CML333 HP301 Il14H Ki3 Ki11 Ky21 M37W M162W Mo18W Ms71 NC350 NC358 Oh7B Oh43 P39 Tx303 Tzi8 EP1 DK105 F7 PE0075", " ")[[1]]
gts31_ph207 = str_split("Mo17 W22 PH207 B97 CML52 CML69 CML103 CML228 CML247 CML277 CML322 CML333 HP301 Il14H Ki3 Ki11 Ky21 M37W M162W Mo18W Ms71 NC350 NC358 Oh7B Oh43 P39 Tx303 Tzi8 EP1 DK105 F7 PE0075", " ")[[1]]
fi = glue("{dirg}2/wheat26.txt")
gts_wheat26 = read_tsv(fi, col_names='gt') %>% pull(gt)

infer_intron <- function(td, zero=T) { # sorted by chrom and start
    #{{{
    chrom = td$chrom[1]; srd = td$srd[1]
    obegs = td$start; oends = td$end
    begs = oends[-nrow(td)]; ends = obegs[-1]
    if(zero == F) begs = begs + 1
    tibble(chrom=chrom, start=begs, end=ends, srd=srd)
    #}}}
}

