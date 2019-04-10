source('functions.R')
gcfg = read_genome_conf()
tg = gcfg$loc.gene %>%
    group_by(gid) %>%
    summarise(chrom = chrom[1], start = min(start), end = max(end)) %>% ungroup() %>%
    mutate(start = as.integer(start), end = as.integer(end))
# read gene mapping
fm = file.path(dirg, 'B73/gene_mapping/maize.v3TOv4.geneIDhistory.txt')
tm = read_tsv(fm, col_names = c("gid", "ngid", "note", "method", "type"))
tm = tm %>% filter(type == '1-to-1') %>% select(ogid=gid,ngid)

read_eqtl_hs <- function(qtag, dird='~/projects/genomes/data') {
    #{{{
    dirw = file.path(dird, qtag)
    fi = sprintf("%s/%s_eQTL.txt", dirw, qtag)
    ti = read_tsv(fi)
    f_hs = sprintf("%s/%s_hs.txt", dirw, qtag)
    t_hs = read_tsv(f_hs)
    # tq: i, gid, qid,qchrom,qstart,qend,qpos, score
    # t_hs: hid,hchrom,hstart,hend, n.tgt
    if(qtag == 'li2013') {
        #{{{
        tqs = ti %>% distinct(qchr.p,qpos.p) %>%
            mutate(qid=sprintf("qtl%04d",1:length(qchr.p)))
        tq = ti %>% select(i=eQTL_ID,gid=`e-traits`,qchrom=qchr.p,qpos=qpos.p,
                           r2=R2) %>%
                inner_join(tqs, by=c('qchrom'="qchr.p",'qpos'="qpos.p")) %>%
                mutate(qstart=qpos, qend=qpos, score=r2) %>%
                select(i, gid, qid,qchrom,qstart,qend,qpos, score)
        t_hs = t_hs %>% mutate(hid=sprintf("hs%03d",1:nrow(t_hs))) %>%
            select(hid,hchrom=Chr,hstart=Start,hend=End,n.tgt=n_gene)
        #}}}
    } else if(qtag == 'liu2017') {
        #{{{
        tqs = ti %>% distinct(eQTL) %>%
            mutate(qid = sprintf("qtl%04d", 1:length(eQTL)))
        tq = ti %>% mutate(i=1:nrow(ti)) %>%
            select(i,gid=Trait,eQTL,r2=Joint_Rsq) %>%
            mutate(locstr = eQTL) %>%
            separate(locstr, c('qchrom','locstr'), sep=':') %>%
            separate(locstr, c('qstart','qend'), sep='\\.\\.') %>%
            mutate(qstart=as.integer(qstart), qend=as.integer(qend),
                   qpos=(qstart+qend)/2) %>%
            inner_join(tqs, by='eQTL') %>% rename(score=r2) %>%
            select(i, gid, qid,qchrom,qstart,qend,qpos, score)
        t_hs = t_hs %>% mutate(hid=sprintf("hs%03d",1:nrow(t_hs))) %>%
            select(hid, hchrom=Chr,hstart=Start,hend=End,n.tgt=N_tgts)
        #}}}
    } else if (qtag == 'wang2018') {
        #{{{
        ti = ti %>% mutate(eQTL=sprintf("%d:%g",eQTL_chr,eQTL_pos))
        tqs = ti %>% distinct(eQTL) %>%
            mutate(qid = sprintf("qtl%04d", 1:length(eQTL)))
        tq = ti %>% mutate(i=1:nrow(ti)) %>%
            select(i,gid=etrait,eQTL,qchrom=eQTL_chr,
                   qstart=left.cM,qend=right.cM,qpos=eQTL_pos,r2=qtl.r2) %>%
            inner_join(tqs, by='eQTL') %>% mutate(score = r2) %>%
            select(i, gid, qid,qchrom,qstart,qend,qpos, score)
        t_hs = t_hs %>% mutate(hid=sprintf("hs%03d",1:nrow(t_hs))) %>%
            select(hid,hchrom=Chr,hstart=start,hend=end,n.tgt=n_genes) %>%
            mutate(hchrom = as.numeric(str_sub(hchrom,4)))
        #{{{ convert genetic position to physical position
        fsz = '~/data/genome/B73/v22/15.sizes'
        tsz = read_tsv(fsz, col_names=c('chrom','psize')) %>%
            filter(chrom %in% sprintf("%d",1:10)) %>%
            mutate(chrom=as.numeric(chrom))
        tq1 = tq %>% group_by(qchrom) %>% summarise(gsize=max(qend)) %>% ungroup()
        tq1 = tq1 %>% inner_join(tsz, by=c('qchrom'='chrom')) %>%
            mutate(rate=psize/gsize) %>% select(qchrom,rate)
        tq = tq %>% inner_join(tq1, by='qchrom') %>%
            mutate(qstart = as.integer(qstart*rate),
                   qend = as.integer(qend*rate),
                   qpos = as.integer(qpos*rate)) %>%
            select(-rate)
        tq = tq %>% mutate(qstart=qpos-1, qend=qpos)
        #}}}
        #}}}
    } else {
        stop('unknown qtag:', qtag)
    }
    list(tq = tq, t_hs = t_hs)
    #}}}
}

qtag = 'li2013'
qtag = 'wang2018'
qtag = 'liu2017'
chrs = sprintf("%d", 1:10)
dirw = file.path(dird, qtag)
res = read_eqtl_hs(qtag)
tq = res$tq; t_hs = res$t_hs

#{{{ # check how cis/trans is defined
te = ti %>% mutate(samchr = gchr.p == qchr.p) %>%
    #filter(type == 'cis') %>%
    mutate(dis = abs(gpos.p - qpos.p))
summary(te$dis)
summary(te %>% filter(dis >= 10000000) %>% pull(dis))

te2 = te %>%
    mutate(type2 = ifelse(gchr.p == qchr.p,
                ifelse(dis <= 10000000, 'cis',
                ifelse(abs(gpos.g-qpos.g) <= 10, 'cis', 'trans')),'trans'))
te2 %>% count(type2, type)
te2 %>% filter(type != type2) %>% select(dis, gpos.g, qpos.g, gchr.g, qchr.g)
#}}}

#{{{ add hotspot info
to1 = t_hs %>% distinct(hchrom,hstart,hend,hid,n.tgt)
to2 = tq %>% distinct(qchrom,qstart,qend,qid)
fo1=file.path(dirw, '01.hs.bed')
write_tsv(to1, fo1, col_names=F)
fo2=file.path(dirw, '01.eqtl.bed')
write_tsv(to2, fo2, col_names=F)
# intersectBed -wao -a 01.hs.bed -b 01.eqtl.bed > 01.hs_eqtl.bed

fi = file.path(dirw, '01.hs_eqtl.bed')
ti = read_tsv(fi, col_names=c('hchrom','hstart','hend','hid','n.tgt',
                              'qchrom','qstart','qend','qid','bp')) %>%
    select(qid,hid,hchrom,hstart,hend,n.tgt) %>%
    filter(qid != '.')
ti %>% count(hid)
#
tq2 = tq %>% left_join(ti, by='qid') %>%
    mutate(hid = ifelse(is.na(hid), qid, hid),
           hchrom=ifelse(is.na(hchrom), qchrom, hchrom),
           hstart=ifelse(is.na(hstart), qstart, hstart),
           hend = ifelse(is.na(hend), qend, hend)) %>%
    select(-qid,-qchrom,-qstart,-qend,-qpos,-n.tgt) %>%
    rename(qid=hid,qchrom=hchrom,qstart=hstart,qend=hend) %>%
    distinct(i,gid,qid,qchrom,qstart,qend,score) %>%
    mutate(qpos = (qstart+qend)/2)
#}}}

#{{{ convert eQTL coordinates from v2 to v4
tb = tq2 %>% inner_join(tm, by = c('gid'='ogid')) %>%
    select(-gid) %>% rename(gid = ngid)

to = tb %>% mutate(qstart=qstart-1) %>%
    distinct(qchrom,qstart,qend,qid) %>%
    arrange(qchrom, qstart)
fo = file.path(dirw, '03.v2.bed')
write_tsv(to, fo, col_names = F)

# liftOver 03.v2.bed chain 04.v4.raw.bed -multiple
# liftOver 04.v4.raw.bed chain 05.v4.bed -multiple

fl = file.path(dirw, "05.v4.bed")
tl = read_tsv(fl, col_names = c("nchr","nbeg","nend",'qid','npiece'), col_types = 'ciici')
tl = tl %>% mutate(size=nend-nbeg) %>%
    arrange(qid, desc(size)) %>%
    group_by(qid) %>%
    summarise(nchr=nchr[1],nbeg=nbeg[1],nend=nend[1]) %>% ungroup()

tb3 = tb %>% inner_join(tl, by = 'qid') %>%
    select(-qchrom,-qstart,-qend,-qpos) %>%
    rename(qchrom=nchr, qstart=nbeg, qend=nend) %>%
    mutate(qpos=(qstart+qend)/2) %>%
    inner_join(tg, by='gid') %>%
    rename(gchrom=chrom, gstart=start, gend=end) %>%
    mutate(gpos=(gstart+gend)/2) %>%
    select(qid,qchrom,qstart,qend,qpos,
        gid,gchrom,gstart,gend,gpos) %>%
    mutate(type = ifelse(qchrom==gchrom & abs(qpos-gpos) < 1000000, 'cis', 'trans'))

hs = tb3 %>% filter(type=='trans') %>% count(qid,qchrom,qstart,qend,qpos) %>%
    rename(n.tgt=n) %>% filter(n.tgt >= 20)
hs.tgt = tb3 %>% filter(type == 'trans', qid %in% hs$qid)

fo = file.path(dirw, '10.rds')
res = list(hs=hs, hs.tgt=hs.tgt)
saveRDS(res, file=fo)
#}}}



#{{{ # visualize
f_size = '~/data/genome/b73/15_intervals/01.chrom.sizes'
t_size = read_tsv(fz, col_names=c('chrom','size'), col_types='ci')
chrs = sprintf("B%02d", 1:10)
tz = flattern_gcoord_prepare(size)

gpos = flattern_gcoord(ti %>% select(chrom=gchrom, pos=gbeg), tz)
qpos = flattern_gcoord(ti %>% select(chrom=qchrom, pos=qpos), tz)
tp = ti %>% mutate(gpos=!!gpos, qpos=!!qpos)

p1 = ggplot(tp) +
  geom_point(aes(x = qpos, y = gpos, color = type), size=0.2) +
  geom_vline(xintercept = tz$beg, alpha=0.3) +
  geom_vline(xintercept = tz$end, alpha=0.3) +
  geom_hline(yintercept = tz$beg, alpha=0.3) +
  geom_hline(yintercept = tz$end, alpha=0.3) +
  scale_x_continuous(name='eQTL position', breaks=tz$pos, labels=tz$chrom, expand=c(0,0)) +
  scale_y_continuous(name='eGene position', breaks=tz$pos, labels=tz$chrom, expand=c(0,0)) +
  scale_color_brewer(palette = "Set1") +
  otheme(xtitle=T, ytitle=T, xtext=T, ytext=T,
         legend.pos='top.center.out', legend.dir='h')
fp = sprintf("%s/11.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 8)
#}}}

#{{{ # identify hotspots - sliding window approach
fi = file.path(dirw, "10.eQTL.v4.tsv")
ti = read_tsv(fi) %>%
    select(qchrom=chrom.v4, qpos=pos.v4, type, gid=gid.v4, gchrom, gbeg, gend) %>%
    mutate(type2 = ifelse(gchrom==qchrom & abs((gbeg+gend)/2-qpos) < 10000000, 'cis', 'trans'))

tq = ti %>%
    transmute(chr=qchrom, beg=qpos-1, end=qpos, id=1:nrow(ti), gid=gid, type=type2) %>%
    arrange(chr, beg)

fo = file.path(dirw, '22.eqtl.bed')
write_tsv(tq, fo, col_names=F)

#{{{
fo = file.path(dirw, '21.gene.bed')
to = tg %>% select(gchrom,gbeg,gend,gid.v4)
write_tsv(to, fo, col_names=F)
#}}}

# winslide.pl -i $genome/Zmays_v4/15.sizes -o 21.win.bed -step 1000000 -size 1000000
# intersectBed -wao -a 21.win.bed -b 22.eqtl.bed > 23.ovlp.bed
# intersectBed -wao -a 21.win.bed -b 21.gene.bed > 21.win.gene.bed

fv = file.path(dirw, "23.ovlp.bed")
tv = read.table(fv, sep="\t", header=F, as.is=T)[,1:9]
colnames(tv) = c("chr","beg","end",'qchr','qbeg','qend','id','gid','type')
tv = as_tibble(tv)

types = c("trans","cis")
tv2 = tv %>% group_by(chr, beg, end) %>%
    summarise(cis = sum(type=='cis'),
              trans=sum(type=='trans'),
              nsnp = length(unique(qend[type=='trans']))) %>% ungroup()
gbeg = flattern_gcoord(tv2 %>% transmute(chrom=chr, pos=beg), t_size)
gend = flattern_gcoord(tv2 %>% transmute(chrom=chr, pos=end), t_size)
tv3 = tv2 %>% mutate(beg=gbeg+1, end=gend, pos=(gbeg+gend)/2+1) %>%
    select(-nsnp) %>%
    gather(type, ngene, -chr, -pos, -beg, -end) %>%
    replace_na(list(ngene=0)) %>%
    mutate(type = factor(type, levels = types))

tp = tv3
p1 = ggplot(tp) +
    #geom_bar(aes(x=pos, y=ngene, fill=type), stat='identity') +
    geom_point(aes(x=pos, y=ngene, color=type), size=.2) +
    geom_segment(data=tz, aes(x=beg,xend=end,y=0,yend=0), size=1.5) +
    #geom_vline(xintercept = tz$beg, alpha=0.3) +
    #geom_vline(xintercept = tz$end, alpha=0.3) +
    scale_x_continuous(name='eQTL position', breaks=tz$pos, labels=tz$chrom, expand=c(0,0)) +
    scale_y_continuous(name='# of eGenes', expand=expand_scale(mult=c(0,.05))) +
    scale_fill_manual(values = brewer.pal(3,"Set1")[2:1]) +
    scale_color_npg() +
    facet_wrap(~type, nrow=2) +
    otheme(xtitle=T,ytitle=T,xtext=T,ytext=T)
fp = sprintf("%s/23.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 6)
#}}}

#{{{ #
fg = file.path(dirw, '21.win.gene.bed')
tg = read.table(fg, sep="\t", header=F, as.is=T)
colnames(tg) = c("chr","beg","end","gchr","gbeg","gend","gid",'bp')
#tg = tg[tg$chr %in% chrs,]
#tg$chr = as.integer(tg$chr)

grp = dplyr::group_by(tg, chr, beg, end)
tg2 = dplyr::summarise(grp, ngene = sum(gid != '.'))
stopifnot(nrow(tg2)==nrow(tv2))

tv5 = cbind.data.frame(tv2, ngene = tg2$ngene)
tv5 = within(tv5, {den = trans/(1*ngene)})
tv5$den[is.na(tv5$den)] = 0

tv6 = tv5[tv5$den >= 1.25 & tv5$trans >= 10,]

fo = file.path(dirw, "25.hs.raw.bed")
write.table(tv6[,1:3], fo, sep = "\t", row.names = F, col.names = F, quote = F)

# mergeBed -i 25.hs.raw.bed > 26.hs.merged.bed
# intersectBed -wao -a 26.hs.merged.bed -b 25.hs.raw.bed > 27.hs.bed

fm = file.path(dirw, '27.hs.bed')
tm = read.table(fm, sep="\t", header=F, as.is=T)[,1:6]
colnames(tm) = c("mchr","mbeg","mend","chr","beg","end")
tms = unique(tm[,1:3])
tms = tms[order(tms$mchr, tms$mbeg),]
tms = cbind(tms, mid = sprintf("hs%02d", 1:nrow(tms)))
tm = merge(tms, tm, by=c("mchr","mbeg","mend"))


tm2 = merge(tm, tv, by=c("chr","beg","end"))[,-c(1:3)]
colnames(tm2)[1:3] = c("chr","beg","end")
tm2 = tm2[order(tm2$chr,tm2$beg),]
fo = file.path(dirw, "29.hs.tsv")
write.table(tm2, fo, sep = "\t", row.names = F, col.names = T, quote = F)

tm3 = tm2[tm2$type=='trans',]
tm4 = data.frame(mid = sprintf("trans-eQTL.%s", tm3$mid), gid = tm3$gid, funcat='', opt='trans-eQTL', stringsAsFactors = F)

fo = file.path(dirw, "29.hs.module.tsv")
write.table(tm4, fo, sep = "\t", row.names = F, col.names = T, na = '', quote = F)

# generate stats
fh = file.path(dirw, "29.hs.tsv")
th = read_tsv(fh)
th %>% group_by(mid) %>% summarise(ntrans=sum(type=='trans')) %>% print(n=70)
#}}}

#{{{ ## identify trans-eQTL hotspots - single SNP approach
fi = file.path(dirw, "10.eQTL.v4.tsv")
ti = read.table(fi, sep="\t", header=T, as.is=T)
ti = within(ti, {type2 = ifelse(gchr==qchr & abs((gbeg+gend)/2-qpos) < 10000000, 'cis', 'trans')})

tq = unique(data.frame(chr=ti$qchr, beg=ti$qpos-1, end=ti$qpos))
tq = tq[order(tq$chr,tq$beg),]

fo = file.path(dirw, '41.eqtl.bed')
write.table(tq, fo, sep = "\t", row.names = F, col.names = F, quote = F)

# bcftools view -H -R 41.eqtl.bed $misc2/mo17vnt/53.vnt.final/02.vcf.gz > 42.vnt.vcf

ti2 = ti[ti$type2 == 'trans',]
ti2 = cbind(ti, direction = sign(ti$ADD))
grp = dplyr::group_by(ti2, qchr, qpos, direction)
ti3 = dplyr::summarise(grp, ntgt = n())
sum(ti3$ntgt >= 20)
ti4 = ti3[ti3$ntgt >= 20,]

ti4 = ti4[ti4$qchr %in% chrs,]
ti4$qchr = as.integer(ti4$qchr)
ti4 = ti4[order(ti4$qchr, ti4$qpos),]
ti4a = ti4[ti4$direction == -1,]
ti4a = cbind(mid = sprintf("trans-eQTL-B.hs%02d", 1:nrow(ti4a)), ti4a)
ti4b = ti4[ti4$direction == 1,]
ti4b = cbind(mid = sprintf("trans-eQTL-M.hs%02d", 1:nrow(ti4b)), ti4b)
ti4 = rbind(ti4a, ti4b)


ti5 = merge(ti2, ti4[,-5], by=c("qchr","qpos",'direction'))
ti5$qchr = as.integer(ti5$qchr)
ti5 = ti5[order(ti5$mid, ti5$gid),]

fo = file.path(dirw, "49.hs.tsv")
write.table(ti5, fo, sep = "\t", row.names = F, col.names = T, na = '', quote = F)

tm = data.frame(mid = ti5$mid, gid = ti5$gid, funcat='', opt='trans-eQTL', stringsAsFactors = F)
fo = file.path(dirw, "49.hs.module.tsv")
write.table(tm, fo, sep = "\t", row.names = F, col.names = T, na = '', quote = F)

# stats
fh = file.path(dirw, "49.hs.module.tsv")
th = read.table(fh, sep="\t", header=T, as.is=T)
ths = ddply(th, .(mid), summarise, size=length(mid))
describe(ths$size)
#}}}



