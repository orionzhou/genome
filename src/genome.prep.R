#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(argparse))

parser <- ArgumentParser(description = 'save genome stats to .rds')
parser$add_argument("genome", default='B73', help="genome")
parser$add_argument("--dirg", default='~/projects/genome/data',
                    help="genome directory [default: %(default)s]")
args <- parser$parse_args()

genome = args$genome
dirg = args$dirg

source("~/projects/genome/src/functions.R")

dirw = file.path(dirg, genome)
fi = file.path(dirw, '15_intervals', '01.chrom.bed')
bed.chrom = read_tsv(fi, col_names = c("chrom", "start", "end")) %>%
    mutate(start = start + 1, end = as.numeric(end))
#
fi = file.path(dirw, '15_intervals', '11.gap.bed')
bed.gap = read_tsv(fi, col_names = c("chrom", "start", "end")) %>%
    mutate(start = start + 1, end = as.numeric(end))
#
fi = file.path(dirw, '15_intervals', '20.gap.sep.60win.tsv')
#bed.60 = read_tsv(fi, col_names = T)
#
fi = file.path(dirw, '50_annotation', '10.tsv')
ti = read_tsv(fi, col_types = 'ccccciic')
t_gs = ti %>%
    filter(etype == 'exon') %>%
    group_by(gid, tid) %>%
    summarise(size = sum(end - start + 1)) %>%
    group_by(gid) %>%
    summarise(size = max(size))
#
fi = file.path(dirw, '50_annotation', '15.tsv')
ti = read_tsv(fi, col_types = 'ccccciic')
loc.gene = ti
#
size.gene = ti %>%
    filter(etype == 'exon') %>%
    group_by(gid, tid) %>%
    summarise(size = sum(end - start + 1)) %>%
    group_by(gid) %>%
    summarise(size = max(size))
cat(sprintf("%5d genes for %s\n", nrow(size.gene), genome))

fi = file.path(dirw, '50_annotation', '15.desc.tsv')
gene.desc = read_tsv(fi)

res = list(chrom=bed.chrom, gap=bed.gap,
           loc.gene=loc.gene, size.gene=size.gene, t_gs=t_gs,
           gene.desc=gene.desc)
fo = file.path(dirw, '55.rds')
saveRDS(res, file = fo)
#e1 = env(bed.chrom = bed.chrom, bed.gap = bed.gap, loc.gene = loc.gene, size.gene = size.gene)
#save(e1, file = fo)

envname = 'genomedata'
#x = attach(file.path(dirg, genome, '55.rda'), name = envname)
#detach(envname, character.only = T)
