source('functions.R')
genome = 'Zmays_B73'
dirw = file.path(dirp, 'data2', genome, 'gene_mapping')
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

fi = file.path(dirw, 'B73_Mo17.tsv')
ti2 = read_tsv(fi, col_names=F) %>% select(gid=X4, Mo17=X8, type=X9)

#{{{ jcvi synmap pipeline output
read_synmap <- function(qry, tgt='B73', diri='/home/springer/zhoux379/projects/wgc/data/raw') {
    #{{{
fi = sprintf("%s/%s_%s/20_synteny/07.t.ortholog", diri, qry, tgt)
read_tsv(fi, col_names=c('tid1','tid2')) %>%
    filter(tid2 != '.') %>%
    mutate(type = ifelse(str_sub(tid2, -1, -1)=="'", 'rbh','syn')) %>%
    mutate(tid2 = str_replace(tid2, "'", '')) %>%
    separate(tid1, c('gid1','iso1'), sep="[\\.\\_]", remove=F) %>%
    separate(tid2, c('gid2','iso2'), sep="[\\.\\_]", remove=F) %>%
    select(gid1, gid2, type, tid1, tid2)
    #}}}
}

qrys = c("Mo17","W22",'PH207')
to = tibble(qry=qrys, tgt='B73') %>%
    mutate(xref = map2(qry, tgt, read_synmap)) %>%
    unnest(xref)

fo = file.path(dirw, 'xref.synmap.tsv')
write_tsv(to, fo)
#}}}


tx %>% count(m1=='', m2=='')
tx %>% filter(m1!='', m2!='') %>% mutate(x=m1==m2) %>% count(x)

tx = ti %>% full_join(ti3, by='gid') %>% rename(m1=Mo17,m2=gid2) %>%
    replace_na(list(m1='',m2=''))
tx %>% count(m1=='', m2=='')
tx %>% filter(m1!='', m2!='') %>% mutate(x=m1==m2) %>% count(x)
