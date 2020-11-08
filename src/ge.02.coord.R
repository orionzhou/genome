require(GenomicRanges)
source('functions.R')
genome = 'Zmays_B73'
gcfg = read_genome_conf(genome)
tsyn = read_syn(gcfg)

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

#{{{ TSS, promoter & regulatory regions
dirw = file.path(dirg, genome, '50_annotation')
fi = file.path(dirw, "15.tsv")
ti = read_tsv(fi) %>% filter(etype == 'exon') %>% mutate(start=start-1)

bp_promoter = 2000
tt = ti %>% filter(ttype=='mRNA') %>%
    group_by(gid, tid) %>%
    summarise(chrom=chrom[1], start=min(start), end=max(end), srd=srd[1]) %>%
    ungroup() %>%
    arrange(chrom, start, end) %>%
    mutate(tss = ifelse(srd=='-', end, start)) %>%
    mutate(pstart = ifelse(srd=='-', tss, tss-bp_promoter)) %>%
    mutate(pend = ifelse(srd=='-', tss+bp_promoter, tss)) %>%
    select(gid, tid, chrom,start,end, srd, tss, pstart, pend)

tt1 = tt %>% mutate(tstart=tss, tend=tss+1, score='.') %>%
    select(chrom,tstart,tend,gid,score,srd) %>% arrange(chrom, tstart)
fo = file.path(dirw, '16.tss.bed')
write_tsv(tt1, fo, col_names = F)

tp = tt %>% mutate(score='.') %>% select(chrom, pstart, pend, gid,score,srd)
fo = file.path(dirw, '16.promoter.2kb.bed')
write_tsv(tt2, fo, col_names = F)
# bedtools getfasta -fi ../10_genome.fna -bed 16.promoter.2kb.bed -name -fo 16.promoter.2kb.fas
# faSize -detailed 16.promoter.2kb.fas > 16.promoter.2kb.sizes
tp = tt %>% mutate(gstart=0,gend=pend-pstart,srd='+', pid=1:n()) %>%
    select(gid,gstart,gend,srd, chrom, pstart, pend, pid)
fo = file.path(dirw, '16.promoter.2kb.chain.bed')
write_tsv(tp, fo, col_names = F)
# chain.py fromBed 16.promoter.2kb.chain.bed 16.promoter.2kb.sizes ../15_intervals/01.chrom.sizes > 16.promoter.2kb.chain

bp_prox = 2000; bp_dist = 1e5
to1 = tt %>% select(chrom,start,end,gid,srd) %>% mutate(opt='genic')
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

fo = file.path(dirw, '16.regulation.bed')
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

#{{{ gene feature intervals for overlap with 100-bp bins
genome = 'Zmays_Mo17'
genome = 'Zmays_PH207'
genome = 'Zmays_W22'
dirw = file.path(dird, genome, '50_annotation')
fi = file.path(dirw, "15.tsv")
#{{{ add introns & tile sub-features
ti = read_tsv(fi)
ti1 = ti %>% filter(ttype == 'mRNA')
ti2 = ti %>% filter(ttype != 'mRNA') %>%
    distinct(chrom,start,end,ttype,etype)
gr = with(ti2, GRanges(seqnames=chrom, ranges=IRanges(start, end=end)))
sum(width(gr)) - sum(width(GenomicRanges::reduce(gr)))
#
tim = ti1 %>% filter(etype=='exon') %>% group_by(gid) %>%
    summarise(chrom=chrom[1], start=min(start), end=max(end), srd=srd[1]) %>%
    ungroup()
tix = ti1 %>% filter(etype=='exon')
tic = ti1 %>% filter(etype=='CDS')
tiu5 = ti1 %>% filter(etype=='five_prime_UTR')
tiu3 = ti1 %>% filter(etype=='three_prime_UTR')
grm = with(tim, GRanges(seqnames=chrom, ranges=IRanges(start, end=end)))
grx = with(tix, GRanges(seqnames=chrom, ranges=IRanges(start, end=end)))
grc = with(tic, GRanges(seqnames=chrom, ranges=IRanges(start, end=end)))
gru5 = with(tiu5, GRanges(seqnames=chrom, ranges=IRanges(start, end=end)))
gru3 = with(tiu3, GRanges(seqnames=chrom, ranges=IRanges(start, end=end)))
grm=GenomicRanges::reduce(grm)
grx=GenomicRanges::reduce(grx)
grc=GenomicRanges::reduce(grc)
gru5=GenomicRanges::reduce(gru5)
gru3=GenomicRanges::reduce(gru3)
setdiff(grx, c(grc,gru5,gru3))
gru5 = setdiff(gru5, grc)
gru3 = setdiff(gru3, c(grc, gru5))
gri = setdiff(grm, c(grc,gru5,gru3))
mcols(grc)['etype'] = 'cds'
mcols(gru5)['etype'] = 'utr5'
mcols(gru3)['etype'] = 'utr3'
mcols(gri)['etype'] = 'intron'
#
x = c(grc,gru5,gru3,gri)
sum(width(x)) - sum(width(GenomicRanges::reduce(x)))
tc = tibble(chrom = as.character(seqnames(x)),
            start = start(x), end = end(x), ttype='mRNA', etype=x$etype) %>%
    bind_rows(ti2) %>% arrange(chrom,start,end) %>%
    mutate(start = start - 1)
#}}}

fo = file.path(dirw, '15.tile.bed')
write_tsv(tc, fo, col_names = F)
#}}}

#{{{ gene bins for meta plots
genome = 'Zmays_B73'
dirw = file.path(dirp, 'data2', genome)
gcfg = read_genome_conf(genome)

tg = gcfg$gene.loc %>% filter(etype=='exon') %>%
    group_by(gid,tid,ttype) %>%
    summarise(chrom=chrom[1],start=min(start)-1,end=max(end),srd=srd[1]) %>%
    ungroup() %>% mutate(ftype='Gene') %>%
    select(fid=gid,ftype,chrom,start,end,srd)
tg %>% count(ftype)

tp = tg %>%
    mutate(len = end - start) %>%
    mutate(start0 = start-len/2, end0 = end+len/2) %>%
    arrange(chrom,start,end) %>%
    mutate(idx=1:n()) %>% filter(len >= 200)
#
gr = with(tp, GRanges(seqnames=chrom, ranges=IRanges(start0,end0)))
nbin = 200
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
    select(chrom, start, end, fid, i, srd) %>%
    arrange(chrom, start, end)

fo = file.path(dirw, '35.metaplot.bin.bed')
write_tsv(tp, fo, col_names=F)
#}}}

#{{{ create bioconductor TxDb
require(GenomicFeatures)
org = 'Zmays_B73'

chromInfo = gcfg$chrom %>% select(chrom,length=size)
f_gff = file.path(dird, org, "50_annotation/10.gff")
txdb = makeTxDbFromGFF(f_gff, format='gff3',
                       organism='Zea mays', chrominfo=chromInfo)

f_txdb = file.path(dird, org, "50_annotation/10.sqlite")
saveDb(txdb, file=f_txdb)

#x = select(txdb, keys=keys(txdb), columns="TXTYPE", keytype="GENEID")
#}}}

