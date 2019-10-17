source('functions.R')
dirw = file.path(dirp, 'data2/TE_annotation')

#{{{ fix Mo17 TE annotation
fm = file.path(dird, 'Zmays_Mo17/download/ids.txt')
tm = read_tsv(fm, col_names='x') %>%
    separate(x, c('oid','x'), sep=' ', extra='merge') %>%
    separate(x, c('x','suf'), sep=',') %>%
    separate(x, c('pre','nid'), sep=' Mo17 ') %>%
    mutate(oid = str_replace(oid, '^>', '')) %>%
    mutate(nid = str_replace(nid, '^chromosome ', '')) %>%
    select(oid, nid)

fi = file.path(dirw, '01.Mo17.gff.gz')
ti = read_tsv(fi, skip=3, col_names=F, col_types='ccciicccc')
to = ti %>% left_join(tm, by=c('X1'='nid')) %>%
    mutate(X1 = oid) %>% select(-oid)
fo = file.path(dirw, '02.Mo17.gff')
write_tsv(to, fo, col_names=F)
#}}}

#{{{ convert TE gff to bed
genome = 'B73'
genome = 'Mo17'
genome = 'PH207'
genome = 'W22'
fi = sprintf("%s/03.%s.gff", dirw, genome)
skip = ifelse(genome %in% c("Mo17",'W22'), 0, 3)
ti = read_tsv(fi, col_names=F, skip=skip)
#
tt = ti %>%
    separate(X9, c('id0','name0'), sep=';') %>%
    separate(id0, c('id0','id'), sep='=') %>%
    separate(name0, c('name0','name'), sep='=') %>%
    select(chrom=X1,src=X2,type=X3,start=X4,end=X5,srd=X7,id,name) %>%
    mutate(sfam = str_sub(id, 0, 3),
        fam = str_sub(id, 4, 8),
        assembly = str_sub(id, 9, 16),
        copy = str_sub(id, 17)) %>%
    select(chrom,start,end,srd,type,sfam,fam,assembly,copy,id,name) %>%
    arrange(chrom,start,end)
nrow(tt) - length(unique(tt$id))
tt %>% count(type,sfam)
tt %>% count(type,sfam,fam)
tt %>% count(assembly)

fo = sprintf("%s/10.%s.tsv", dirw, genome)
write_tsv(tt, fo)
tb = tt %>% mutate(start=start-1)
fo = sprintf("%s/11.%s.bed", dirw, genome)
write_tsv(tt, fo, col_names=F)
#}}}


