require(GenomicRanges)
source('functions.R')
genome = 'Zmays_B73'
gcfg = read_genome_conf(genome)
tsyn = read_syn(gcfg)

#{{{ exon intervals for each gene
dirw = file.path(dird, genome, '50_annotation')
fi = file.path(dirw, "10.tsv")
ti = read_tsv(fi) %>% filter(etype == 'exon')

gr = with(ti, GRanges(seqnames=chrom, ranges=IRanges(start, end=end), gid=gid))
grl = reduce(split(gr, elementMetadata(gr)$gid))
x = unlist(grl)
tc = tibble(gid = names(x), chrom = as.character(seqnames(x)),
            start = start(x), end = end(x))

tc %>% group_by(gid) %>%
    summarise(size = sum(end - start + 1)) %>%
    ungroup() %>% group_by(1) %>%
    summarise(n = n(), me = mean(size), md = median(size),
              min = min(size), max = max(size),
              q25 = quantile(size, .25), q75 = quantile(size, .75))

ta = tc %>% mutate(start = start - 1) %>% select(chrom, start, end, gid) %>%
    arrange(chrom, start)
fo = file.path(dirw, '10.ase.bed')
write_tsv(ta, fo, col_names = F)
#}}}

#{{{ CDS intervals for pop-gene analysis
dirw = file.path(dird, genome, '50_annotation')
fi = file.path(dirw, "15.tsv")
ti = read_tsv(fi) %>% filter(etype == 'CDS')

to = ti %>% mutate(loc = sprintf("%s:%d-%d", chrom, start, end)) %>%
    distinct(loc)
fo = file.path(dirw, '15.cds.txt')
write_tsv(to, fo, col_names = F)

tid = ti %>% distinct(gid)
fo = file.path(dirw, '15.gid.txt')
write_tsv(tid, fo, col_names = F)
#}}}

#{{{ TSS for each gene
dirw = file.path(dirg, genome, '50_annotation')
fi = file.path(dirw, "15.tsv")
ti = read_tsv(fi) %>% filter(etype == 'exon') %>% mutate(start=start-1)

tt = ti %>% filter(ttype=='mRNA') %>%
    group_by(gid, tid) %>%
    summarise(chrom=chrom[1], start=min(start), end=max(end), srd=srd[1]) %>%
    ungroup() %>%
    mutate(tss = ifelse(srd=='-', end, start)) %>%
    arrange(chrom, tss) %>%
    mutate(gstart=start, gend=end) %>%
    mutate(start=tss, end=tss+1) %>%
    select(chrom,start,end,srd, gid, tid, gstart,gend)

fo = file.path(dirw, '15.tss.bed')
write_tsv(tt, fo, col_names = F)

bp_prox = 2000; bp_dist = 1e5
to1 = tt %>% select(chrom,start=gstart,end=gend,gid,srd) %>% mutate(opt='genic')
to2 = tt %>%
    mutate(s1 = ifelse(srd=='-', start, start-bp_prox)) %>%
    mutate(e1 = ifelse(srd=='-', start+bp_prox, start)) %>%
    select(chrom,start=s1,end=e1,srd,gid) %>% mutate(opt='proximal')
to3 = tt %>%
    mutate(s1 = ifelse(srd=='-', start+bp_prox, start-bp_dist)) %>%
    mutate(e1 = ifelse(srd=='-', start+bp_dist, start-bp_prox)) %>%
    select(chrom,start=s1,end=e1,srd,gid) %>% mutate(opt='distal')
to = rbind(to1, to2, to3) %>% arrange(chrom,start,end) %>%
    filter(end > 0) %>% mutate(start = pmax(start,0))

#fo = file.path(dirw, '15.promoter.bed')
#write_tsv(to2, fo, col_names = F)
fo = file.path(dirw, '15.regulation.bed')
write_tsv(to, fo, col_names = F)
#}}}

#{{{ gap stats
genomes = c("B73", "W22", "PH207")

tp = tibble()
for (org in orgs) {
	fi = file.path(dirg, org, "15_intervals/11.gap.bed")
	ti = read_tsv(fi, col_names = c("chr", "beg", "end"))
	tp = rbind(tp, ti %>% mutate(org = !!org) %>% select(org, everything))
}
tp = tp %>% mutate(size = end - beg)

tps = tp %>%
    group_by(org) %>%
    summarise(grp, num = n(), me = mean(size), md = median(size),
              min = min(size), max = max(size), total = sum(size))
#}}}

x = read_tsv(file.path(dirw, 'x.bed'), col_names=F)
x %>% group_by(X5) %>% summarise(bp = sum(X10)) %>% ungroup() %>%
    count(bp >= 200)

#{{{ gene and TE intervals for bs-seq analysis
genome = 'Zmays_B73'
dirw = file.path(dird, genome, '50_annotation')
gcfg = read_genome_conf(genome)
tg = gcfg$gene.loc %>% filter(etype=='exon') %>%
    mutate(ftype='Gene') %>% mutate(start=start-1) %>%
    select(fid=gid, ftype, etype, chrom, start, end, srd)
#
ft = file.path(dirw, '30.TE.tsv')
tt = read_tsv(ft) %>% mutate(ftype='TE') %>%
    select(fid=id,ftype,chrom,start,end,srd) %>%
    mutate(start = start-1, etype='exon') %>%
    select(fid, ftype, etype, everything())
#
tx = rbind(tg, tt)

#{{{ get up, down, TE-exon intervals
txs = tx %>%
    group_by(fid, ftype) %>%
    summarise(chrom=chrom[1],start=min(start),end=max(end),srd=srd[1]) %>%
    ungroup() %>% mutate(off = floor((end-start)/2))
#
txu = txs %>%
    mutate(start1 = ifelse(srd=='-', end, start-off),
           end1 = ifelse(srd=='-', end+off, start), etype='up') %>%
    mutate(eid = sprintf("%s|%s|%s", ftype,fid,etype)) %>%
    select(eid,chrom,start=start1,end=end1,srd)
txd = txs %>%
    mutate(start1 = ifelse(srd=='-', start-off, end),
           end1 = ifelse(srd=='-', start, end+off), etype='down') %>%
    mutate(eid = sprintf("%s|%s|%s", ftype,fid,etype)) %>%
    select(eid,chrom,start=start1,end=end1,srd)
txt = tx %>% filter(ftype == 'TE') %>%
    mutate(eid = sprintf("%s|%s|%s", ftype,fid,etype)) %>%
    select(eid,chrom,start,end,srd)

to = rbind(txu, txd, txt) %>% filter(start >= 0, end >= 0) %>%
    filter(end-start >= 10) %>% mutate(srd = ifelse(srd=='.', '*', srd))
gr = with(to, GRanges(seqnames=chrom, ranges=IRanges(start+1,end), strand=srd))
grl = tile(gr, n=10)
names(grl) = to$eid

tvi = as_tibble(grl) %>%
    select(eid=group_name,chrom=seqnames,start,end,srd=strand) %>%
    mutate(chrom=as.character(chrom))
tv1 = tvi %>% filter(srd != '-') %>%
    arrange(eid, start) %>%
    group_by(eid) %>%
    mutate(i = 1:length(eid)) %>% ungroup()
tv2 = tvi %>% filter(srd == '-') %>%
    arrange(eid, desc(start)) %>%
    group_by(eid) %>%
    mutate(i = 1:length(eid)) %>% ungroup()
tva = rbind(tv1, tv2) %>% mutate(start = start-1)
#}}}

#{{{ create exon/intron intervals
# get intron intervals (only for genes with >=2 exons)
gids_2exon = tg %>% count(fid) %>% filter(n>1) %>% pull(fid)
txi = tg %>% filter(fid %in% gids_2exon) %>%
    arrange(chrom, start) %>%
    group_by(fid) %>% nest() %>%
    mutate(intron=map(data, infer_intron)) %>%
    select(-data) %>% unnest() %>% mutate(ftype='Gene',etype='intron') %>%
    select(fid,ftype,etype,chrom,start,end,srd)

tp = rbind(tg, txi) %>%
    mutate(eid = sprintf("%s|%s|%s", ftype,fid,etype)) %>%
    select(eid,chrom,start,end,srd) %>%
    mutate(size = end-start) %>%
    filter(size >= 10) %>% mutate(srd = ifelse(srd=='.', '*', srd)) %>%
    arrange(eid, start)

to = tp %>%
    group_by(eid, chrom, srd) %>%
    summarise(chromStart=start[1], chromEnd=end[length(eid)],
        thickStart=chromStart, thickEnd=chromEnd,
        blockCount = n(), blockSizes = str_c(size, collapse=','),
        blockStarts = str_c(start-chromStart, collapse=',')) %>%
    ungroup() %>% mutate(score=0,itemRgb=0) %>%
    select(chrom,chromStart,chromEnd,eid,score,srd,thickStart,thickEnd,
        itemRgb,blockCount, blockSizes, blockStarts)

fo = file.path(dirw, 't1.bed')
write_tsv(tp, fo, col_names=F)

tp1 = tp %>% filter(srd != '-') %>%
    arrange(eid, chrom, start) %>%
    group_by(eid) %>%
    mutate(i = 1:length(eid)) %>% ungroup()
tp2 = tp %>% filter(srd == '-') %>%
    arrange(eid, chrom, desc(start)) %>%
    group_by(eid) %>%
    mutate(i = 1:length(eid)) %>% ungroup()
to1 = rbind(tp1, tp2)

# make chain bed
to2 = to1 %>% mutate(len=end-start) %>%
    arrange(eid,i) %>%
    group_by(eid) %>%
    mutate(rstart = c(0,cumsum(len)[-length(eid)])) %>%
    mutate(rend = rstart+len) %>% ungroup() %>%
    select(eid, rstart, rend, srd, chrom, start, end) %>%
    mutate(idx = 1:length(eid))
tos = to2 %>% group_by(eid) %>% summarise(size=max(rend), srd=srd[1]) %>%
    ungroup() %>% mutate(srd=ifelse(srd=='.','*',srd))

fo = file.path(dirw, 't1.bed')
write_tsv(to2, fo, col_names=F)
fo = file.path(dirw, 't2.sizes')
write_tsv(tos, fo, col_name = F)
# chain.py fromBed t1.bed t2.sizes ../15_intervals/01.chrom.sizes > tmp.chain

# tile exon/intron intervals
tps = tp %>% group_by(eid) %>% summarise(size=sum(size), srd=srd[1]) %>%
    ungroup() %>% mutate(srd=ifelse(srd=='.','*',srd)) %>%
    filter(size >= 10)
gr = with(tps, GRanges(seqnames=eid, ranges=IRanges(1, size), strand=srd))
grl = tile(gr, n=10)
names(grl) = tps$eid

tvi = as_tibble(grl) %>% select(eid=seqnames,start,end,srd=strand)
tvo = tvi %>% mutate(start=start-1) %>%
    arrange(eid, start) %>%
    group_by(eid) %>%
    mutate(i = 1:length(eid)) %>% ungroup()

to = tvo %>% mutate(chrom=eid, score='.') %>%
    mutate(eid=sprintf("%s|%s",eid,i)) %>%
    select(chrom,start,end,eid,score,srd)
fo = file.path(dirw, 'tmp.i.bed')
write_tsv(to, fo, col_names=F)
# liftOver -multiple -bedPlus=6 tmp.i.bed tmp.chain tmp.o.bed tmp.unmap

#}}}

fi = file.path(dirw, 'tmp.o.bed')
ti = read_tsv(fi, col_names=c('chrom','start','end','eid','eidx','srd'))
tx = ti %>%
    separate(eid, c('type','gid','part','i'), sep='[\\|]') %>%
    unite('eid', type:part, sep="|") %>%
    select(eid, chrom, start, end, srd, i)

tp = rbind(tva, tx) %>% select(chrom,start,end,srd, eid,i) %>%
    mutate(chrom = as.character(chrom)) %>%
    arrange(chrom, start, end)
fo = file.path(dirw, '35.intervals.bed')
write_tsv(tp, fo, col_names=F)
#}}}

#{{{ [obsolete] gene and TE intervals for bs-seq analysis
genome = 'Zmays_B73'
dirw = file.path(dird, genome, '50_annotation')
gcfg = read_genome_conf(genome)

tg = gcfg$gene.loc %>% filter(etype=='exon') %>%
    group_by(gid,tid,ttype) %>%
    summarise(chrom=chrom[1],start=min(start)-1,end=max(end),srd=srd[1]) %>%
    ungroup() %>% mutate(ftype='Gene') %>%
    select(fid=gid,ftype,chrom,start,end,srd)
tg %>% count(ftype)
#
ft = file.path(dirw, '30.TE.tsv')
tt = read_tsv(ft) %>% mutate(ftype='TE') %>%
    select(fid=id,ftype,chrom,start,end,srd)

tp = rbind(tg, tt) %>%
    mutate(len = end - start) %>%
    mutate(start0 = start-len/2, end0 = end+len/2) %>%
    arrange(chrom,start,end) %>%
    group_by(1) %>% mutate(idx=1:n()) %>% ungroup() %>% select(-`1`)
#
gr = with(tp, GRanges(seqnames=chrom, ranges=IRanges(start0,end0)))
nbin = 40
grl = tile(gr, n=nbin)

tp1 = tp %>% select(idx,fid,ftype,srd)
to = as_tibble(grl) %>% select(idx=group,chrom=seqnames,start,end) %>%
    inner_join(tp1, by='idx') %>%
    select(-idx)
#
to1 = to %>% filter(srd != '-') %>%
    arrange(fid, chrom, start) %>%
    group_by(fid) %>%
    mutate(i = 1:n()) %>% ungroup()
to2 = to %>% filter(srd == '-') %>%
    arrange(fid, chrom, desc(start)) %>%
    group_by(fid) %>%
    mutate(i = 1:n()) %>% ungroup()
tp = rbind(to1, to2) %>%
    mutate(start = start-1) %>%
    filter(start>=0, end>=0) %>%
    select(chrom, start, end, fid, ftype, i) %>%
    arrange(chrom, start, end)

fo = file.path(dirw, '35.intervals.bed')
write_tsv(tp, fo, col_names=F)
#}}}

#{{{ ACR
# liftOver -bedPlus=3 acr.0.bed ../08_seq_map/mapf.chain acr.1.bed unmap
# sortBed -i acr.1.bed | cut -f1-3 > acr.2.bed
# intersectBed -a ../50_annotation/15.regulation.bed -b acr.2.bed -wao > acr.3.bed
## closestBed -a acr.2.bed -b ../50_annotation/15.tss.bed -d > acr.3.bed
dirw = file.path(dird, genome, 'chromatin')
opts = c('promoter','gene_body','distal','none')
bp_acr = 50
fa = file.path(dirw, 'acr.3.bed')
ta = read_tsv(fa, col_names=c("chrom",'start','end','gid','srd','opt',
                              'chrom2','start2','end2','bp')) %>%
    group_by(gid,opt) %>% summarise(bp = sum(bp)) %>% ungroup() %>%
    mutate(ovlp = bp >= bp_acr) %>%
    select(gid, opt, ovlp) %>% spread(opt, ovlp) %>%
    mutate(opt = ifelse(distal, 'distal', 'none')) %>%
    mutate(opt = ifelse(genic, 'gene_body', opt)) %>%
    mutate(opt = ifelse(proximal, 'promoter', opt)) %>%
    mutate(opt = factor(opt, levels=opts)) %>%
    select(gid, opt)
ta %>% count(opt)

ta %>% inner_join(tsyn, by='gid') %>%
    group_by(ftype, opt) %>%
    summarise(ng = n()) %>%
    mutate(prop = ng/sum(ng)) %>%
    ungroup() %>% select(ftype, opt, prop) %>% spread(opt, prop)
#}}}

#{{{ UMR
# awk -F'\t' '$4 == "Unmethylated"' umr.1.bed | cut -f1-3 > umr.2.bed
# intersectBed -a ../50_annotation/15.regulation.bed -b umr.2.bed -wao > umr.3.bed
opts = c('promoter','gene_body','distal','none')
bp_umr = 200
fu = file.path(dirw, 'umr.3.bed')
tu = read_tsv(fu, col_names=c("chrom",'start','end','gid','srd','opt',
                              'chrom2','start2','end2','bp')) %>%
    group_by(gid,opt) %>% summarise(bp = sum(bp)) %>% ungroup() %>%
    mutate(ovlp = bp >= bp_umr) %>%
    select(gid, opt, ovlp) %>% spread(opt, ovlp) %>%
    mutate(opt = ifelse(distal, 'distal', 'none')) %>%
    mutate(opt = ifelse(genic, 'gene_body', opt)) %>%
    mutate(opt = ifelse(proximal, 'promoter', opt)) %>%
    mutate(opt = factor(opt, levels=opts)) %>%
    select(gid, opt)
tu %>% count(opt)

tu %>% inner_join(tsyn, by='gid') %>%
    group_by(ftype, opt) %>%
    summarise(ng = n()) %>%
    mutate(prop = ng/sum(ng)) %>%
    ungroup() %>% select(ftype, opt, prop) %>% spread(opt, prop)

to = ta %>% rename(acr = opt) %>% inner_join(tu, by='gid') %>% rename(umr=opt)
to %>% count(acr, umr)
fo = file.path(dirw, 'chromatin.tsv')
write_tsv(to, fo)
#}}}



