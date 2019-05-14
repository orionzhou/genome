#{{{
require(devtools)
load_all('~/git/rmaize')
require(tidyverse)
dirp = '~/projects/genome'
dird = file.path(dirp, 'data2')
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
fi = '~/projects/genome/data/B73/50_annotation/15.tsv'
ti = read_tsv(fi)
#
infer_intron <- function(td) {
    #{{{
    chrom = td$chrom[1]; srd = td$srd[1]
    begs = td$start; ends = td$end
    tibble(chrom = chrom,
           start = ends[-nrow(td)]+1, end=begs[-1],
           srd = srd)
    #}}}
}

tx = ti %>% filter(etype=='exon') %>% select(-etype)
gids_2exon = tx %>% count(gid) %>% filter(n>1) %>% pull(gid)
t_intr = tx %>% filter(gid %in% gids_2exon) %>%
    arrange(chrom, start) %>%
    group_by(gid) %>% nest() %>%
    mutate(intron=map(data, infer_intron)) %>%
    select(-data) %>% unnest() %>% mutate(etype = 'intron')

etype_map = c("CDS"='cds','five_prime_UTR'='utr5','three_prime_UTR'='utr3')
to = ti %>% filter(ttype == 'mRNA', etype != 'exon') %>%
    mutate(etype=etype_map[etype]) %>%
    select(gid, etype, chrom, start, end, srd)
to = rbind(to, t_intr)
to1 = to %>% filter(srd == '+') %>%
    arrange(chrom, start) %>%
    group_by(gid, etype) %>%
    mutate(eidx = 1:n()) %>% ungroup()
to2 = to %>% filter(srd == '-') %>%
    arrange(chrom, desc(start)) %>%
    group_by(gid, etype) %>%
    mutate(eidx = 1:n()) %>% ungroup()

to = rbind(to1, to2) %>% arrange(chrom, start) %>%
    mutate(start = start - 1) %>%
    select(chrom,start,end,srd,gid,etype,eidx)
fo = file.path(dirw, "12.genic.bed")
write_tsv(to, fo, col_names = F)

# intersectBed -a 11 -b 12 -wao > 13

etypes = c('cds','utr5','utr3','intron')
fi = file.path(dirw, "13.ovlp.bed")
ti = read_tsv(fi, col_names = c('mchrom','mbeg','mend','mid','sids',
                                'echrom','ebeg','eend','esrd','gid',
                                'etype','eidx','bp')) %>%
    select(mid, sids, gid, etype, eidx, bp) %>%
    mutate(etype = factor(etype, levels=etypes)) %>%
    arrange(mid, sids, gid, etype, desc(bp), eidx) %>%
    group_by(mid, sids, gid) %>%
    summarise(etype=etype[1], eidx=eidx[1]) %>% ungroup()
ti %>% mutate(genic = ifelse(gid == '.', T, F)) %>% count(genic)
mu = ti %>% filter(gid != '.')

fo = file.path(dirw, '15.mu.genic.tsv')
write_tsv(mu, fo)
#}}}

#{{{ add meta
fi = file.path(dirw, "15.mu.genic.tsv")
mu = read_tsv(fi)
n_mus = c('1','2','>=3')
tg1 = mu %>%
    #filter(etype != 'intron') %>%
    filter(etype == 'cds') %>%
    group_by(gid) %>%
    summarise(n_mu = n(),
              mid = str_c(mid, sep = "|", collapse = "|"),
              sids = str_c(sids, sep = "|", collapse = "|")) %>% ungroup() %>%
    mutate(n_mu_p = ifelse(n_mu >= 3, ">=3", n_mu)) %>%
    mutate(n_mu_p = factor(n_mu_p, levels = n_mus))

# add TF info
ft = file.path(dirp, 'data/B73/61_functional/06.tf.tsv')
tt = read_tsv(ft)
tts = tt %>% count(fam) %>% rename(fam_size=n)
tt = tt %>% inner_join(tts, by='fam') %>%
    group_by(fam) %>% mutate(fam_idx = 1:n()) %>% ungroup() %>%
    mutate(fam_idx_size = sprintf("%d/%d", fam_idx, fam_size)) %>%
    mutate(tf = T)
tg2 = tg1 %>% left_join(tt, by = 'gid') %>%
    replace_na(list(tf = F))

# add W22 gene ID
fi = '~/projects/genome/data/B73/gene_mapping/B73_W22.tsv'
ti = read_tsv(fi,col_names=c('chr1','start1','end1','gid1','chr2','start2','end2','gid2','src','type')) %>%
    transmute(gid=gid1, gid_W22=gid2, map_type=type)
ti2 = ti %>% group_by(gid) %>%
    summarise(gid_W22 = paste(gid_W22,collapse=','),
              map_type = paste(sort(unique(map_type)), collapse=',')) %>% ungroup()
tg3 = tg2 %>% left_join(ti2, by='gid')

# add W22-B73 gene model change
fi = '~/projects/wgc/data/05_stats/10.B73_W22.tsv'
ti = read_tsv(fi)
impacts = c("no_change","low",'modifier','moderate','high','non-syntenic')
tg4 = tg3 %>% left_join(ti, by = 'gid') %>%
    replace_na(list(syn='non-syntenic')) %>%
    mutate(impact = ifelse(syn=='syntenic', impact, syn)) %>%
    mutate(impact = factor(impact, levels = impacts)) %>%
    select(-po, -syn,-tid)
tg4 %>% filter(tf) %>% count(map_type, impact) %>% print(n=30)

# output
to = tg4
fo = file.path(dirw, "16.gene.mu.tsv")
write_tsv(to, fo)

to %>% count(n_mu_p)
to %>% count(impact)
to %>% count(tf)
low_impacts = c("no_change","low",'modifier','moderate')
to %>% filter(impact %in% low_impacts)
to %>% filter(tf, n_mu != '1', impact %in% low_impacts)
to %>% filter(tf, n_mu != '1', impact %in% low_impacts)
#}}}

#{{{ process W22 expression
fi1 = file.path(dirw, 'Samples_B73_Ref_HTseq.txt')
fi2 = file.path(dirw, 'Samples_W22_Ref_HTseq.txt')
ti = read_tsv(fi1) %>% gather(cond, counts, -Genes) %>%
    rename(gid=Genes) %>%
    filter(!str_detect(gid, '^_'))
#
ti2 = ti %>% mutate(cond = str_replace(cond, "_Ref_counts.txt", "")) %>%
    separate(cond, c('cond','ref'), sep='_') %>%
    separate(cond, c('Genotype',"Tissue",'Rep'), sep='-') %>%
    group_by(gid,Genotype,Tissue) %>%
    summarize(counts = sum(counts)) %>% ungroup() %>%
    group_by(Genotype, Tissue) %>%
    mutate(rpm = counts/sum(counts) * 1000000) %>%
    ungroup()
#
tis_map = c('A'='Anther','En'='Endosperm','Em'='Embryo','I'='Internode',
'IE'='Ear', 'L'='Leaf','L10'='Leaf10', 'R'='Root', 'SC'='Shoot', 'T'='Tassel')
ti3 = ti2 %>% filter(Genotype %in% c("B","W")) %>%
    filter(Tissue %in% c("En","Em","L10","SC","R","I")) %>%
    mutate(Tissue=tis_map[Tissue]) %>%
    mutate(cond = sprintf("%s_%s", Genotype, Tissue)) %>%
    select(gid,cond,rpm) %>%
    mutate(rpm = sprintf("%.01f", rpm)) %>%
    spread(cond, rpm)
te = ti3
#}}}

#{{{ further characterize all TF mutants
fi = file.path(dirw, "16.gene.mu.tsv")
ta = read_tsv(fi)

# add TF45 info
fi = '~/projects/grn/data/08_y1h/10.tsv'
ti = read_tsv(fi) %>% distinct(reg) %>% transmute(gid=reg, TF45=T)
ta2 = ta %>% left_join(ti, by='gid') %>% replace_na(list(TF45=F))
ta2 = ta2 %>% filter(tf | TF45) %>% select(-tf)

# add biomap support
fi = '~/projects/grn/data/14_eval_sum/02.valid.bm.spc.rds'
ti = readRDS(fi)$tf
ta3 = ta2 %>% left_join(ti, by=c('gid'='reg.gid')) %>%
	replace_na(list(n.tgt=0))

# add eQTL support
fi = '~/projects/grn/data/14_eval_sum/02.hs.tsv'
ti = read_tsv(fi)
ta4 = ta3 %>% left_join(ti, by=c('gid'='reg.gid')) %>%
	rename(eQTL=qtags, eQTL_n=n_qtag, eQTL_fc=fc, eQTL_grp_size=max.grp.size) %>%
	replace_na(list(eQTL=''))

# add W22 expression
ta5 = ta4 %>% left_join(te, by='gid') %>%
    mutate(max_W22_exp=pmax(W_Embryo,W_Endosperm,W_Internode,W_Leaf10,W_Root,W_Shoot))

to = ta5
to = to %>% rename(gid_B73=gid) %>%
    select(-n_mu_p, -fam_idx_size, -ttype)
fo = file.path(dirw, '20.tf.tsv')
write_tsv(to, fo)
#}}}

#{{{ select TF mutants
fi = file.path(dirw, '20.tf.tsv')
ti = read_tsv(fi)
ti
ti2 = ti %>% filter(n_mu >= 1, map_type == 'One-to-One', max_W22_exp >= 2)

gids1 = ti2 %>% filter(eQTL != '') %>%
    #filter(str_detect(eQTL,',')) %>%
    pull(gid_B73)
gids2 = ti2 %>% filter(n.tgt >= 3) %>% arrange(desc(n.tgt)) %>%
    filter(row_number() <= 20) %>% pull(gid_B73)
gids3 = ti2 %>%
    filter(fam %in% c("HSF","LBD","SBP","TCP","WRKY","MYB"), n_mu >= 2) %>%
    pull(gid_B73)
to1 = tibble(select='eQTL', gid=gids1)
to2 = tibble(select='biomAP', gid=gids2)
to3 = tibble(select='multi-fam', gid=gids3)
to = rbind(to1,to2,to3) %>% group_by(gid) %>%
    summarise(select=paste(select,collapse=',')) %>% ungroup()
to = ti %>% inner_join(to, by=c('gid_B73'='gid')) %>%
    select(gid_B73, select, everything()) %>%
    arrange(select, fam, gid_B73)

fo = file.path(dirw, '30.tf.selected.tsv')
write_tsv(to, fo)
#}}}

#{{{ stock order
fi = file.path(dirw, '30.tf.selected.tsv')
ti = read_tsv(fi)
fi = file.path(dirw, "15.mu.genic.tsv")
mu = read_tsv(fi)
etypes = c("cds")
fi = file.path(dirw, "31.erika.selected.xlsx")
tk = read_xlsx(fi) %>% select(gid=gid_B73, mu) %>%
    replace_na(list(mu='')) %>%
    mutate(mu = map(mu, str_split, "[\\|]")) %>% mutate(mu = map(mu, 1)) %>%
    unnest()

tm = ti %>% select(gid=gid_B73, select_reason=select, n_mu=n_mu) %>%
    inner_join(mu, by='gid') %>% filter(etype %in% etypes)
tm %>% count(gid, n_mu) %>% mutate(nd = n-n_mu) %>% pull(nd)
# use Erika's list to filter
tm.1 = tm %>% filter(!gid %in% tk$gid)
tm.2 = tm %>% inner_join(tk, by=c('gid'='gid','mid'='mu'))
tm2 = rbind(tm.1, tm.2)
#
tm2 = tm2 %>% mutate(etype=fct_relevel(etype, etypes)) %>%
    arrange(gid, etype, eidx) %>%
    group_by(gid) %>%
    filter(row_number() <= 3) %>% ungroup()
tm3 = tm2 %>% select(-n_mu) %>%
    mutate(sid = str_split(sids, ',')) %>%
    mutate(sid = map_chr(sid, 1)) %>%
    arrange(select_reason, gid)

fo = file.path(dirw, '32.gene.stocks.tsv')
write_tsv(tm3, fo)

tmo = tm3 %>% group_by(sid) %>%
    summarise(mid=paste(mid, collapse=","),
              gid=paste(gid, collapse=','),
              select_reason = paste(select_reason, collapse=',')) %>%
    ungroup() %>%
    arrange(sid)

fo = file.path(dirw, '34.stocks.tsv')
write_tsv(tmo, fo)
#}}}

#{{{ gather info from existing/collaborator stocks
fi = file.path(dirw, '34.stocks.tsv')
ti = read_tsv(fi)

# previous stocks
f_od = file.path(dirw, 'Springer_UniformMu_orders.xlsx')
t_od = read_xlsx(f_od, col_names = c("row",'stock','season')) %>%
    mutate(note = sprintf("%s_%s",season, row)) %>%
    group_by(stock) %>%
    summarise(note = paste(note, collapse=",")) %>%
    ungroup()

# mgc reply
fi = file.path(dirw, '40.stocks.misc.xlsx')
t_mgc = read_xlsx(fi, col_names = c("sid",'note')) %>%
    mutate(note = factor(note, unique(note))) %>%
    arrange(sid, desc(note)) %>%
    group_by(sid) %>% summarise(note = note[1]) %>% ungroup()

to = ti %>%
    left_join(t_mgc, by='sid') %>% rename(note1=note) %>%
    left_join(t_od, by=c('sid'='stock')) %>% rename(note2=note) %>%
    arrange(note1,note2,sid)
to %>% count(is.na(note1), is.na(note2))

fo = file.path(dirw, '41.stocks.field.tsv')
write_tsv(to, fo, na='')
#}}}
