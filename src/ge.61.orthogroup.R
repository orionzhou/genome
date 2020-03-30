source('functions.R')
dirw = file.path(dirp, 'data2/ortholog')
gcfg = read_genome_conf()

fi = file.path(dirw, '01.maize.arath.tsv')
ti = read_tsv(fi, col_names=c("gid1",'gid2','code','genus','fam')) %>%
    separate(gid1,c('org1','gid1','uid'), sep='[|]') %>%
    separate(gid1, c('src1','gid1'), sep='=') %>%
    separate(gid2,c('org2','gid2','uid'), sep='[|]') %>%
    separate(gid2, c('src2','gid2'), sep='=')

ti2 = ti %>% filter(gid1 %in% gcfg$gene$gid)
ti2 %>% count(org1, src1, org2, src2)
ti3 = ti2 %>% select(gid1, gid2, code, fam) %>%
    arrange(gid1, gid2, code) %>%
    group_by(gid1, gid2) %>% slice(1) %>% ungroup()

gids1 = ti3 %>% count(gid1) %>% filter(n>1) %>% pull(gid1)
gids2 = ti3 %>% count(gid2) %>% filter(n>1) %>% pull(gid2)
ti4 = ti3 %>% mutate(tag1=gid1 %in% gids1, tag2=gid2 %in% gids2) %>%
    mutate(tag=ifelse(tag1, ifelse(tag2, 'mm', '1m'), ifelse(tag2, 'm1', '11'))) %>%
    select(-tag1,-tag2)
ti4 %>% count(tag)

to1 = ti4 %>% filter(tag=='11')
to2 = ti4 %>% filter(tag=='1m') %>% arrange(gid1,code) %>% group_by(gid1) %>% slice(1) %>% ungroup()
to3 = ti4 %>% filter(tag=='m1') %>% arrange(gid2,code) %>% group_by(gid2) %>% slice(1) %>% ungroup()
to4 = ti4 %>% filter(tag=='mm') %>% filter(code=='LDO')

to = rbind(to1,to2,to3,to4)
to %>% count(code, tag)

fo = file.path(dirw, '10.maize.arath.tsv')
write_tsv(to, fo)
