Maize UniformMu mutants characterization
================
November 02, 2018





































Raw GFF3 file mapping to B73 AGPv3 were obtained from maizeGDB - 79,407
insertion sites and related stocks information

Insertion sites were lifted over to AGPv4 - 77,085 sites

The overlap of these sites with AGPv4 exons were checked: - 36,648 sites
overlap w. at least one exon(s) - 38,778 sites do not overlap w. an exon
(i.e., either in introns, UTRs or intergenic spaces)

A total of 14,831 genes have at least one exonic insertion sites, 9,046
have at least two exonic insertions.

Table 1. Summary of genes with at least one exonic UniformMu insertion
sites.

| \# exonic insertions | \# genes |
| -------------------: | -------: |
|                    1 |     5784 |
|                    2 |     3347 |
|                 \>=3 |     5699 |

Structural integrity of these 14,831 B73 genes in the W22 genome were
checked based on variants called by whole genome comparison between B73
and W22. 10,236 genes have either no, low (such as synonymous change) or
moderate (such as mis-sense) levels of changes in W22 compared to B73.

Table 2. Summary of gene changes between B73 and W22.

|       impact | \# genes |
| -----------: | -------: |
|   no\_change |     2102 |
|          low |     1190 |
|     modifier |     1770 |
|     moderate |     5174 |
|         high |     3247 |
| non-syntenic |     1347 |

The expression of these genes in 23 tissues of B73 and Mo17 were
checked, with genes showing differential expression between B73 and Mo17
in at least 5 tissues prioritized.

Table 3. Summary of genes showing DE in different number of tissues
between B73 and Mo17.

| \# tissues DE btw. B73 and Mo17 | \# genes |
| ------------------------------: | -------: |
|                               0 |     1775 |
|                             1-5 |     5097 |
|                            5-10 |     4798 |
|                             10+ |     3160 |

Finally, 826 of these genes are transcription factor based on
[PlantTFDB](http://planttfdb.cbi.pku.edu.cn/index.php?sp=Zma).

[This table](/data/uniformmu/15.uniformmu.exon.tsv) lists all 14,831
genes with at least one exonic UniformMu insertion sites, with columns:
- `gid`: AGPv4 gene ID - `n_mu`: number exonic UniformMu insertion sites
- `mid`, `sids`: mutant ID and stock ID - `tf`: whether this gene is TF
- `ttype`: type of gene, ncRNA or protein-coding mRNA - `impact`, `eff`:
impact and effect of the change from B73 to W22 - `de_BM`: number of
tissues showing DE between B73 and Mo17

A total of 339 TFs have at least 2 exonic insertions and have gene
structure conserved in the W22 genome.

Table 4. Summary of the 339 TFs by number DE tissues between B73 and
Mo17

| \# tissues DE btw. B73 and Mo17 | \# genes |
| ------------------------------: | -------: |
|                               0 |       32 |
|                             1-5 |      134 |
|                            5-10 |      122 |
|                             10+ |       51 |
