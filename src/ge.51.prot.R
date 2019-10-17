require(GenomicRanges)
source('functions.R')
genome = 'Zmays_B73'
gcfg = read_genome_conf(genome)

dird = file.path(dirp, 'data2')
dirw = file.path(dird, '51_proteome')

# bcftools consensus -s sample -f $ref/10_genome.fna $re/data/ase_bcf/MO17.bcf > 10.genome.Mo17.fna
gt = "Mo17"
