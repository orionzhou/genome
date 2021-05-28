source('functions.R')
dirw = glue('{dirp}/data2/syntelog')
gcfg = read_genome_conf()

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


#{{{ jcvi synmap pipeline output
read_synmap2 <- function(qry, tgt='B73', diri='~/projects/wgc/data/raw') {
    #{{{
fi = sprintf("%s/Zmays_%s-Zmays_%s/xref.t.tsv", diri, qry, tgt)
read_tsv(fi, col_names=c('tid1','tid2')) %>%
    filter(tid2 != '.') %>%
    mutate(type = ifelse(str_sub(tid2, -1, -1)=="'", 'rbh','syn')) %>%
    mutate(tid2 = str_replace(tid2, "'", '')) %>%
    separate(tid1, c('gid1','iso1'), sep="[\\.\\_]", remove=F) %>%
    separate(tid2, c('gid2','iso2'), sep="[\\.\\_]", remove=F) %>%
    select(gid1, gid2, type, tid1, tid2)
    #}}}
}
read_synmap <- function(qry, tgt='B73', diri='~/projects/wgc/data/raw') {
    #{{{
    fi = glue("{diri}/Zmays_{qry}-Zmays_{tgt}/xref.t.tsv")
    read_tsv(fi, col_names=c('tid1','tid2', 'type')) %>%
        filter(tid2 != '.') %>%
        separate(tid1, c('gid1','iso1'), sep="[\\.\\_]", remove=F) %>%
        separate(tid2, c('gid2','iso2'), sep="[\\.\\_]", remove=F) %>%
        select(gid1, gid2, type, tid1, tid2)
    #}}}
}

qrys = c("Mo17","W22",'PH207')
qrys = gts31_ph207
to = tibble(qry=qrys, tgt='B73') %>%
    mutate(xref = map2(qry, tgt, read_synmap)) %>%
    unnest(xref)
to %>% count(qry,tgt) %>% print(n=40)

fo = glue('{dirw}/xref.maize.tsv')
write_tsv(to, fo)
#}}}

#{{{ make xref gene model tibble
fi = glue('{dirw}/xref.maize.tsv')
ti = read_tsv(fi)

gts = c("B73",gts31_ph207)
tg = tibble(gt=gts) %>% mutate(fi=glue("{dirg}/Zmays_{gt}/50_annotation/15.tsv")) %>%
    mutate(x = map(fi, read_tsv)) %>%
    select(gt, x) %>% unnest(x)

tg1 = tg %>% filter(gt=='B73') %>% select(-tid)
tg2 = ti %>% inner_join(tg, by=c('qry'='gt','gid2'='gid','tid2'='tid')) %>%
    select(gt=qry,gid=gid1,gid2,type,ttype,etype,chrom,start,end,srd)
to = tg1 %>% bind_rows(tg2)

top = to %>% filter(srd=='+')
tx = top %>% group_by(gt,gid) %>% summarise(tss=min(start)) %>% ungroup()
top = top %>% inner_join(tx, by=c('gt','gid')) %>%
    mutate(beg=start-tss+1, end=end-tss+1)

tom = to %>% filter(srd=='-')
tx = tom %>% group_by(gt,gid) %>% summarise(tss=max(end)) %>% ungroup()
tom = tom %>% inner_join(tx, by=c('gt','gid')) %>%
    mutate(beg=tss-end+1, end=tss-start+1)

to = top %>% bind_rows(tom) %>% select(gt,gid,type=etype,ttype,beg,end,gid2) %>%
    arrange(gt, gid, type, beg)
fo = glue("{dirw}/maize.gene.model.rds")
saveRDS(to, fo)
#}}}




