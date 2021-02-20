source('functions.R')
dirw = glue("{dirg}/raw/Sviridis_ME034v")

fi = glue("{dirw}/raw.gtf")
ti = read_tsv(fi, col_names=F)

ireplace <- function(s) str_replace(str_replace_all(s,".* ", ''), "\"", "")
ti2 = ti %>% separate(X9, c('gid','tid','enum','eid'), sep="; ") %>%
    mutate(gid = str_replace_all(gid, ".* ", '')) %>%
    mutate(gid = str_replace_all(gid, "\"", '')) %>%
    mutate(tid = str_replace_all(tid, ".* ", '')) %>%
    mutate(tid = str_replace_all(tid, "\"", '')) %>%
    mutate(eid = str_replace_all(eid, ".* ", '')) %>%
    mutate(eid = str_replace_all(eid, "\"", '')) %>%
    mutate(gid = str_replace_all(gid, "-.*", '')) %>%
    mutate(tid = str_replace(tid, "-mRNA-", '.')) %>%
    mutate(eid = str_replace(tid, "-mRNA-", '.'))

tig = ti2 %>% group_by(gid) %>%
    summarise(X1=X1[1], X2=X2[1], X3='gene', X4=min(X4), X5=max(X5),
        X6='.', X7=X7[1], X8='.') %>% ungroup()

to = ti2 %>% bind_rows(tig) %>%
    arrange(gid, X4,X5) %>%
    mutate(note = glue('gene_id \"{gid}\"; transcript_id \"{tid}\"; exon_id \"{eid}\"')) %>%
    mutate(note = noquote(note)) %>%
    select(-gid,-tid,-enum,-eid)


fo = glue("{dirw}/raw2.gtf")
write.table(to, file=fo, quote=F, sep="\t", row.names=F, col.names=F)
#write_tsv(to, fo, na='', quote_escape=F, col_names=F)
