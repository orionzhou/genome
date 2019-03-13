source('functions.R')
source(file.path(dirr, 'Location.R'))
dirg = '~/data/genome/Zmays_v4'

#{{{ create region file using large gaps for gatk
dirw = file.path(dirg, "17_regions")
ft = file.path(dirg, "15.bed")
tt = read_tsv(ft, col_names = F, col_types = 'cii') %>%
    transmute(chrom = X1, start = X2+1, end = X3, size = X3-X2)
fp = file.path(dirg, "16.gap.bed")
tp = read_tsv(fp, col_names = F, col_types = 'cii') %>%
    transmute(chrom = X1, start = X2+1, end = X3, size = X3-X2) %>%
    filter(size >= 50000)

grt = with(tt, GRanges(seqnames = chrom, ranges = IRanges(start, end = end)))
grp = with(tp, GRanges(seqnames = chrom, ranges = IRanges(start, end = end)))
grn = setdiff(grt, grp)
tn = as_tibble(grn) %>%
    mutate(chrom = as.character(seqnames)) %>%
    select(chrom, start, end) %>%
    mutate(size = end - start + 1, gap = '')

to1 = tp %>% mutate(gap = 'gap') %>%
    bind_rows(tn) %>%
    filter(chrom %in% sprintf("%d", 1:10)) %>%
    mutate(chrom = as.integer(chrom)) %>%
    arrange(chrom, start)
fo1 = file.path(dirw, "01.tsv")
write_tsv(to1, fo1)

fi = file.path(dirw, "03.hand.curated.tsv")
ti = read_tsv(fi) %>%
    filter(dbreak == 'y')
grp = with(ti, GRanges(seqnames = chrom, ranges = IRanges(start, end = end)))
grn = setdiff(grt, grp)
tn = as_tibble(grn) %>%
    mutate(chrom = as.character(seqnames)) %>%
    select(chrom, start, end) %>%
    mutate(size = end - start + 1) %>%
    filter(chrom %in% sprintf("%d", 1:10)) %>%
    mutate(rid = sprintf("r%02d", 1:length(chrom)))
summary(tn$size)
tn %>% top_n(-4, size)

fn = file.path(dirw, "10.gap.sep.59win.tsv")
write_tsv(tn, fn)
#}}}


