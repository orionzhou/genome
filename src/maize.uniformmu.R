#{{{
require(tidyverse)
dirp = '~/projects/genomes'
dird = file.path(dirp, 'data')
dirw = file.path(dird, 'uniformmu')
#}}}

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
fi = '~/data/genome/B73/50_annotation/10.tsv'
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
to %>% mutate(exonic = ifelse(gid == '.', T, F)) %>% count(exonic)

tg1 = to %>% filter(gid != '.') %>%
    group_by(gid) %>% 
    summarise(n_mu = n(), 
              mid = str_c(mid, sep = "|", collapse = "|"),
              sids = str_c(sids, sep = "|", collapse = "|")) %>% ungroup() %>%
    mutate(n_mu = ifelse(n_mu >= 3, ">=3", n_mu)) 

# add TF info
ft = '~/data/genome/B73/61_functional/06.tf.tsv'
tt = read_tsv(ft) %>% distinct(gid) %>% mutate(tf = T)
tg2 = tg1 %>% left_join(tt, by = 'gid') %>%
    replace_na(list(tf = F))

# add W22-B73 lookup table
fi = '~/projects/wgc/data/05_stats/10.B73_W22.tsv'
ti = read_tsv(fi)
tg3 = tg2 %>% left_join(ti, by = 'gid') %>%
    filter(!is.na(syn)) %>%
    mutate(impact = ifelse(syn=='syntenic', impact, syn)) %>%
    select(-po, -syn,-tid)

# add B73-Mo17 variation data
fi = '~/projects/briggs/data/42_de/11.de.dom.rda'
x = load(fi)
tm0 = tm %>% select(gid, Tissue, pDE) %>% 
    filter(!is.na(pDE) & pDE != 'non_DE') %>%
    count(gid) %>% 
    rename(n_de_BM = n) %>%
    mutate(de_BM = cut(n_de_BM, breaks = c(0,1,5,10,Inf),
                       labels = c("0",'1-5','5-10','10+'), right = F)) %>%
    select(-n_de_BM)
tg4 = tg3 %>% left_join(tm0, by = 'gid') %>% replace_na(list(de_BM='0'))

tp = tg4
fo = file.path(dirw, "15.uniformmu.exon.tsv")
write_tsv(tp, fo)

tp %>% count(n_mu)
impacts = c("no_change","low",'modifier','moderate')
tp %>% filter(impact %in% impacts) %>% count(de_BM)
tp %>% filter(tf, n_mu >= 2, impact %in% impacts) %>% count(de_BM)
#}}}



