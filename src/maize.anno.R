#{{{
require(tidyverse)
#
dirw = '~/data/genome/B73/61_functional'
fg = file.path(dirw, "../50_annotation/10.gtb")
tg = read_tsv(fg, col_types = 'ccciiccccccccccccc')
gids = unique(tg$par)
length(gids)
#
fm = file.path(dirw, "../gene_mapping/maize.v3TOv4.geneIDhistory.txt")
tm = read_tsv(fm, col_names = F) %>%
    transmute(ogid = X1, gid = X2, change = X3, method = X4, type = X5)
length(unique(tm$gid))
sum(unique(tm$gid) %in% gids)
#
Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}
#}}}

#{{{ # interpro
dirg = '~/data/genome/Zmays_v4/61.interpro'

fg = file.path(dirg, "../51.gtb")
tg = read.table(fg, sep = "\t", header = T, as.is = T)[,1:2]
colnames(tg) = c("tid", "gid")

fm = file.path(dirg, "../57.longest.tsv")
tm = read.table(fm, sep = "\t", header = F, as.is = T)
colnames(tm) = c("gid", 'tid')

ft = file.path(dirg, "07.tsv")
tt = read.table(ft, sep = "\t", header = F, as.is = T)
colnames(tt) = c("tid", "gos")

tt2 = merge(tg, tt, by = 'tid', all.x = T)
stopifnot(nrow(tg) == nrow(tt2))
tt2 = tt2[order(tt2$gid, tt2$tid),]

fo1 = file.path(dirg, "10_mrna.tsv")
write.table(tt2[,c(1,3)], fo1, sep = "\t", row.names = F, col.names = F, quote = F, na = '')

tm2 = merge(tm, tt, by = 'tid', all.x = T)
stopifnot(nrow(tm2) == nrow(tm))
tm2 = tm2[order(tm2$gid, tm2$tid),]
fo2 = file.path(dirg, "11_gene.tsv")
write.table(tm2[,c(2,3)], fo2, sep = "\t", row.names = F, col.names = F, quote = F, na = '')
#}}}

#{{{ maize-GAMER V4 GO
fgo = '~/data/db/go/go-basic.tsv'
tgo = read_tsv(fgo) %>%
    transmute(goid = goid, level = level, goname = name)

to = tibble()
ctags = c("Interproscan5", "arabidopsis", "argot2.5", "fanngo", 
          "pannzer", "uniprot.plants", "aggregate")
for (ctag in ctags) {
    fi = sprintf("%s/raw/maize.B73.AGPv4.%s.gaf", dirw, ctag)
    ti = read_tsv(fi, skip = 1)[,c(2,5,7,9)]
    colnames(ti) = c("gid", "goid", "evidence", "gotype")
    ti = ti %>% mutate(ctag = ctag) %>%
        select(ctag, gid, goid, evidence, gotype)
    to = to %>% bind_rows(ti)
}

goids = unique(to$goid)
length(goids)
sum(goids %in% tgo$goid)
to2 = to %>% inner_join(tgo, by = 'goid')
to2 %>% count(ctag)
to2 %>% distinct(ctag, gid) %>% count(ctag)
to2 %>% count(ctag, evidence)

fo = file.path(dirw, "01.go.tsv")
write_tsv(to2, fo)
#}}}

#{{{ maize-GAMER v3 gold-standard GO
fi = sprintf("%s/raw/maize_v3.gold.gaf", dirw)
ti = read_tsv(fi, skip = 1)[,c(3,5,7,9)]
colnames(ti) = c("ogid", "goid", "evidence", "gotype")

to = ti %>%
    inner_join(tm, by = 'ogid') %>%
    filter(gid %in% gids) %>% 
    filter(type == '1-to-1') %>%
    select(goid, gid, evidence, gotype)
length(unique(to$gid))
length(unique(to$goid))

fo = file.path(dirw, '02.go.gs.tsv')
write_tsv(to, fo)
#}}}

#{{{ v3 TFDB to v4
fi = file.path(dirw, 'raw/Zma_TF_list')
ti = read_tsv(fi, col_names = T) %>%
    transmute(ogid = Gene_ID, otid = TF_ID, fam = Family) %>%
    filter(ogid != '') %>%
    group_by(ogid) %>%
    summarise(fam = Mode(fam))

tf = ti %>%
    inner_join(tm, by = 'ogid') %>%
    filter(gid %in% gids) %>%
    #filter(type == '1-to-1') %>%
    select(gid, fam)
tfs = tf %>% distinct(fam) %>% arrange(fam) %>%
    mutate(fid = sprintf("tf%04d", 1:length(fam)))
tf = tf %>% inner_join(tfs, by = 'fam') %>%
    select(fid, gid, fam) %>%
    arrange(fid, gid)

tf %>% count(fid)
fo = file.path(dirw, "06.tf.tsv")
write_tsv(tf, fo)
#}}}

#{{{ # housekeeping genes
dirw = file.path(Sys.getenv("genome"), "Zmays_v4/housekeeping")

fi = file.path(dirw, '01.tsv')
ti = read.table(fi, header = F, sep = "\t", as.is = T)
colnames(ti) = c("ogid", 'name')

ti2 = merge(ti, tm, by = 'ogid')
sum(ti2$ngid %in% tg$par)

fo = file.path(dirw, "11.tsv")
write.table(ti2, fo, sep = "\t", row.names = F, col.names = T, quote = F)
#}}}

#{{{ # CornCyc@maizeGDB query-results (obsolete)
fi = file.path(dirw, "query-results.tsv")
ti = read_tsv(fi, col_names = c('pathway', 'gids'))

ptn = "([&;\"])|(</?((i)|(sub)|(sup))>)"
to = ti %>% 
    mutate(pathway = str_replace_all(pathway, ptn, ""),
           gids = str_replace_all(gids, "[()\"]", "")) %>%
    mutate(gid = str_split(gids, " ")) %>%
    select(-gids) %>%
    mutate(pid = sprintf("cc%04d", 1:length(pathway))) %>%
    unnest() %>%
    select(pid, gid, pathway)
#}}}

#{{{ PMN CornCyc
fi = file.path(dirw, "raw/corncyc_pathways.20180702")
ti = read_tsv(fi)

to = ti %>% transmute(pid = `Pathway-id`,
                      pname = `Pathway-name`,
                      gid = `Gene-name`) %>%
    filter(gid != 'unknown', gid %in% gids) %>%
    distinct(pid, pname, gid) %>%
    select(pid, gid, pname)

fo = file.path(dirw, '07.corncyc.tsv')
write_tsv(to, fo)
#}}}

#{{{ PPIM
fi = file.path(dirw, 'raw/highConfidentPPIs.txt')
ti = read_tsv(fi, col_names = F) %>%
    transmute(tid1 = X1, tid2 = X2)

tt = tibble(tid = unique(c(ti$tid1, ti$tid2))) %>%
    mutate(strtype = str_detect(tid, "GRMZM"))
tt1 = tt %>% filter(strtype) %>% select(-strtype) %>%
    separate(tid, c("gid", "suf"), sep = "_", remove = F) %>%
    select(-suf)
tt2 = tt %>% filter(!strtype) %>% select(-strtype) %>%
    mutate(gid = str_replace(tid, "_FGP", "_FG"))
tt = tt1 %>% bind_rows(tt2)

ti = ti %>% inner_join(tt, by = c("tid1" = "tid")) %>%
    rename(gid1 = gid) %>%
    inner_join(tt, by = c("tid2" = 'tid')) %>%
    rename(gid2 = gid) %>%
    select(-tid1, -tid2)

to = ti %>% rename(ogid1 = gid1, ogid2 = gid2) %>%
    inner_join(tm, by = c("ogid1" = "ogid")) %>%
    rename(type1 = type, gid1 = gid) %>%
    select(ogid1, ogid2, gid1, type1) %>%
    inner_join(tm, by = c("ogid2" = "ogid")) %>%
    rename(type2 = type, gid2 = gid) %>%
    select(gid1, gid2, type1, type2)
to %>% count(type1, type2) %>% print(n=25)

fo = file.path(dirw, '08.ppim.tsv')
write_tsv(to, fo)
#}}}

#{{{ merge functional annotations from different sources
fi = file.path(dirw, "01.tsv")
ti = read_tsv(fi)

fp = file.path(dirw, '../07.corncyc.tsv')
tp = read_tsv(fp) %>%
    transmute(ctag = 'corncyc', gid = gid, goid = pid,
              evidence = '', gotype = '', level = '', goname = pathway)

ff = file.path(dirw, '../TF/10.tsv')
tf = read_tsv(ff) %>%
    transmute(ctag = 'tfdb', gid = gid, goid = fid,
              evidence = '', gotype = '', level = '', goname = fam)

to = rbind(ti, tp, tf)

to %>% count(ctag)
to %>% distinct(ctag, gid) %>% count(ctag)
fo = file.path(dirw, "10.tsv")
write_tsv(to, fo)
#}}}



