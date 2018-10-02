#{{{
require(tidyverse)
dirp = '~/projects/genomes'
dird = file.path(dirp, 'data')
dirw = file.path(dird, 'uniformmu')

#{{{ reformat chr labels
fi = file.path(dirw, "uniformmu_latest_release_v3.gff3")
ti = read_tsv(fi, col_names = F) 
to = ti %>% mutate(X1 = str_replace(X1, "Chr([\\d+])", "\\1"))
fo = file.path(dirw, "01.uniformmu_v3.gff")
write_tsv(to, fo, col_names = F)
#}}}

# liftover v3 -> v4 -> B01 coord

#{{{ create UniformMu Bed
fi = file.path(dirw, "03.v4.gff")
ti = read_tsv(fi, col_names = F) %>%
    mutate(chrom = X1, beg = X4, end = X5, note = X9) %>%
    separate(note, c("namestr", "stockstr"), sep = ";") %>%
    separate(namestr, c("namet", "name"), sep = "=") %>%
    separate(stockstr, c("stockt", "stocks"), sep = "=") %>%
    select(chrom, beg, end, name, stocks) %>%
    arrange(chrom, beg, end)

to = ti %>% mutate(beg = beg - 1) 
fo = file.path(dirw, "11.uniformmu.bed")
write_tsv(to, fo, col_names = F)
#}}}

#{{{ create exon range BED
fi = '/home/springer/zhoux379/data/genome/B73/50_annotation/10.tsv'
ti = read_tsv(fi)

to = ti %>% filter(etype == 'exon') %>%
    select(chrom, start, end, srd, gid, tid) %>%
    mutate(start = start - 1) %>%
    arrange(chrom, start)
fo = file.path(dirw, "12.exon.bed")
write_tsv(to, fo, col_names = F)
#}}}
# intersectBed -a 11 -b 12 -wao > 13

#{{{
fi = file.path(dirw, "13.ovlp.bed")
ti = read_tsv(fi, col_names = F) %>%
    mutate(mchrom = X1, mbeg = X2, mend = X3, mid = X4, sids = X5,
           echrom = X6, ebeg = X7, eend = X8, esrd = X9, gid = X10, tid = X11,
           bp = X12) %>%
    select(mid, sids, gid, tid, bp) 

to = ti %>% group_by(mid, sids, gid) %>%
    summarise(bp = max(bp)) %>% ungroup()

ft = '/home/springer/zhoux379/data/genome/B73/61_functional/06.tf.tsv'
tt = read_tsv(ft) %>% distinct(gid) %>% mutate(tf = T)

tz = to %>% left_join(tt, by = 'gid')
tz %>% dplyr::count(tf)

tz %>% filter(gid != '.') %>%
    group_by(gid) %>% summarise(n_mu = n()) %>% ungroup() %>%
    count(n_mu)

tz %>% filter(gid != '.', !is.na(tf)) %>%
    group_by(gid) %>% summarise(n_mu = n()) %>% ungroup() %>%
    count(n_mu)
#}}}

#}}}
