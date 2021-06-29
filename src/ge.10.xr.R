source('functions.R')
dirw = glue('{dirp}/data2/syntelog')
gcfg = read_genome_conf()
#{{{ functions
read_synmap <- function(qry, tgt='Zmays_B73', diri='~/projects/wgc/data/raw') {
    #{{{
    fi = glue("{diri}/{qry}-{tgt}/xref.t.tsv")
    ti = read_tsv(fi, col_names=c('tid1','tid2', 'type')) %>%
        filter(tid2 != '.')
    if (str_detect(qry, '^Atauschii_')) {
        ti = ti %>% separate(tid2, c('gid2','iso2'), sep="[\\.]", remove=F)
    } else {
        ti = ti %>% separate(tid2, c('gid2','iso2'), sep="[\\.\\_]", remove=F)
    }
    if (tgt == 'Atauschii_AS60') {
        ti = ti %>% separate(tid1, c('iso1','gid1'), sep="[\\.]", remove=F)
    } else {
        ti = ti %>% separate(tid1, c('gid1','iso1'), sep="[\\.\\_]", remove=F)
    }
    ti %>% select(gid1, gid2, type, tid1, tid2)
    #}}}
}
#}}}

#{{{ maizeGDB xref
fi = file.path(dirw, 'gene_model_xref_v4.txt')
ti = read_tsv(fi, skip=4)
xref = ti %>% slice(-nrow(ti)) %>%
    select(gid=1,tid=6,
           sorghum=sorghum_ortholog,
           foxtail_millet=foxtail_millet_ortholog,
           rice=rice_ortholog,
           brachy=brachypodium_ortholog,
           arabidopsis=arabidopsis_ortholog,
           b73=`B73(Zm00001d.2)`,
           W22=`W22(Zm00004b.1)`,
           B104=`B104(Zm00007a.1)`,
           PH207=`PH207(Zm00008a.1)`,
           Mo17a=`Mo17(Zm00009a.1)`,
           EP1=`EP1(Zm00010a.1)`,
           F7=`F7(Zm00011a.1)`,
           Mo17=`Mo17(Zm00014a.1)`)

fo = file.path(dirw, 'xref.maizeGDB.tsv')
write_tsv(xref, fo, na='')
#}}}


#{{{ syntelog xref table
org='maize'; gt0="Zmays_B73"; v='v4'
org='maize'; gt0="Zmays_B73v5"; v='v5'
org='wheat'; gt0="Taestivum_D"; v='csd'
#org='wheat'; gt0="Atauschii_AS60"; v='as60'

qrys = if(org=='maize') glue("Zmays_{gts31_ph207}") else gts_wheat26[!gts_wheat26 %in% gt0]
gts = c(gt0,qrys)
to = tibble(qry=qrys, tgt=gt0) %>%
    mutate(xref = map2(qry, tgt, read_synmap)) %>%
    unnest(xref)
to %>% count(qry,tgt,type) %>% spread(type,n) %>% print(n=40)

fo = glue('{dirw}/xref.{org}.{v}.tsv')
write_tsv(to, fo)
#}}}

#{{{ syntelog xref gene structure tibble
ti = to
fi = glue('{dirw}/xref.{org}.{v}.tsv')
ti = read_tsv(fi)

tg = tibble(gt=gts) %>% mutate(fi=glue("{dirg}/{gt}/50_annotation/15.tsv")) %>%
    mutate(x = map(fi, read_tsv)) %>%
    select(gt, x) %>% unnest(x)
tg2 = tg %>% group_by(gt,gid,tid,ttype) %>% nest() %>% ungroup()
#
tg3a = tg2 %>% filter(gt==gt0) %>% select(-tid, -ttype)
tg3b = ti %>% inner_join(tg2, by=c('qry'='gt','gid2'='gid','tid2'='tid')) %>%
    select(gt=qry,gid=gid1,gid2,type,data)
#
to = tg3a %>% bind_rows(tg3b) %>%
    select(gt, gid, type, gid2, gene=data)
to %>% count(gt,type) %>% spread(type,n) %>% print(n=40)
#
fo = glue("{dirw}/xref.{org}.{v}.rds")
saveRDS(to, fo)
#}}}

#{{{ # obsolete TSS-based locations
tg1 = tg %>% filter(gt==gt0) %>% select(-tid) %>% rename(cstart=start,cend=end)
tg2 = ti %>% inner_join(tg, by=c('qry'='gt','gid2'='gid','tid2'='tid')) %>%
    select(gt=qry,gid=gid1,gid2,type,ttype,etype,chrom,cstart=start,cend=end,srd)
to = tg1 %>% bind_rows(tg2)
#
top = to %>% filter(srd=='+')
tx = top %>% group_by(gt,gid) %>% summarise(tss=min(cstart)) %>% ungroup()
top = top %>% inner_join(tx, by=c('gt','gid')) %>%
    mutate(beg=cstart-tss+1, end=cend-tss+1)
#
tom = to %>% filter(srd=='-')
tx = tom %>% group_by(gt,gid) %>% summarise(tss=max(cend)) %>% ungroup()
tom = tom %>% inner_join(tx, by=c('gt','gid')) %>%
    mutate(beg=tss-cend+1, end=tss-cstart+1)

to = top %>% bind_rows(tom) %>%
    select(gt,gid,type=etype,ttype,beg,end,gid2,chrom,cstart,cend,srd) %>%
    arrange(gt, gid, type, beg)
#}}}


