source('functions.R')
require(GenomicRanges)
genome = 'B73'

#{{{ exon intervals for each gene
dirw = file.path(dirg, genome, '50_annotation')
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


