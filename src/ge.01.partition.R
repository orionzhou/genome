source('functions.R')
require(GenomicRanges)
dirg = file.path(dird, 'Zmays_B73')

dirw = file.path(dirg, "15_intervals")
#{{{ create region file using large gaps for gatk
ft = file.path(dirw, "01.chrom.bed")
tt = read_tsv(ft, col_names = F, col_types = 'cii') %>%
    transmute(chrom = X1, start = X2+1, end = X3, size = X3-X2)
fp = file.path(dirw, "11.gap.bed")
tp = read_tsv(fp, col_names = F, col_types = 'cii') %>%
    transmute(chrom = X1, start = X2+1, end = X3, size = X3-X2)

assign_bin <- function(vs, bin_size=5e7) {
    #{{{
    bins = c()
    cur_size = 0
    bin = 1
    for (v in vs) {
        if(cur_size + v <= bin_size) {
            bins = c(bins, bin)
            cur_size = cur_size + v
        } else {
            if(cur_size > 0) bin = bin + 1
            bins = c(bins, bin)
            cur_size = v
        }
    }
    bins
    #}}}
}

bin_size = 2e7; min_gap_size = 5e3
bin_size = 5e7; min_gap_size = 30e3
bin_size = 1e9; min_gap_size = 5e3
tp1 = tp %>% filter(size >= min_gap_size)
grt = with(tt, GRanges(seqnames = chrom, ranges = IRanges(start, end = end)))
grp = with(tp1, GRanges(seqnames = chrom, ranges = IRanges(start, end = end)))
grn = setdiff(grt, grp)
tn = as_tibble(grn) %>%
    mutate(chrom = as.character(seqnames)) %>%
    select(chrom, start, end) %>%
    mutate(size = end - start + 1) %>%
    group_by(chrom) %>% mutate(bin = assign_bin(size, bin_size=bin_size)) %>%
    ungroup()
tw = tn %>%
    group_by(chrom, bin) %>%
    summarise(start=min(start), end=max(end), size=sum(size)) %>%
    ungroup()
fmt = "r%02d"
if (nrow(tw) >= 100) fmt = "r%03d"
tw = tw %>% mutate(rid = sprintf(fmt, 1:length(chrom))) %>% select(-bin) %>%
    mutate(region=sprintf("%s:%d-%d", chrom, start, end))
nrow(tw)
tw %>% arrange(size) %>% print(n=10)
tw %>% arrange(-size) %>% print(n=10)
tw %>% print(n=10)

fw = sprintf("%s/20.win%d.tsv", dirw, nrow(tw))
write_tsv(tw, fw)
#}}}


