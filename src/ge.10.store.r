source('functions.R')

prepare_genome_data <- function(genome = 'B73', dirg = '~/data/genome') {
    #{{{
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
    fi = file.path(dirw, '50_annotation', '15.tsv')
    ti = read_tsv(fi, col_types = 'ccccciiccc')
    loc.gene = ti
    loc.exon = ti %>% filter(etype == 'exon') %>% select(-etype)
    #
    size.gene = ti %>%
        filter(etype == 'exon') %>%
        group_by(gid, tid) %>%
        summarise(size = sum(end - start + 1)) %>%
        group_by(gid) %>%
        summarise(size = max(size))
    #
    cat(sprintf("%5d genes for %s\n", nrow(size.gene), genome))
    fo = file.path(dirw, '55.rda')
    save(bed.chrom, bed.gap, loc.gene, size.gene, file = fo)
    #e1 = env(bed.chrom = bed.chrom, bed.gap = bed.gap, loc.exon = loc.exon, size.gene = size.gene)
    #save(e1, file = fo)
    T
    #}}}
}

genomes = c("B73", "Mo17", "W22", "PH207", "PHB47")
map_lgl(genomes, prepare_genome_data)

genome = 'B73'
#prepare_genome_data(genome)
envname = 'genomedata'
x = attach(file.path(dirg, genome, '55.rda'), name = envname)
detach(envname, character.only = T)
