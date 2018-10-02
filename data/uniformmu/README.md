# maize UniformMu insertion sites and stock information

Raw GFF3 file mapping to B73 AGPv3 were obtained from maizeGDB
  - 79,407 insertion sites and related stocks information

Insertion sites were lifted over to AGPv4
  - 77,085 sites

The overlap of these sites with AGPv4 exons were checked:
  - 36,648 sites overlap w. at least one exon(s)
  - 38,778 sites do not overlap w. an exon (i.e., either in introns, UTRs or intergenic spaces)

A total of 14,831 genes have at least one exonic insertion sites, 9,046 have at least two exonic insertions.

| # exonic insertions | # genes |
|----------:|-------------:|
| 1 |  5,785 |
| 2 |  3,347 |
| >=3 | 5,699 |

Table containing mu insertion site IDs and affected gene IDs: [link](15.uniformmu.exon.tsv)

Of these affected genes, a total of 826 TFs have at least one exonic insertion sites, 498 have at least two exonic insertions.

| # exonic insertions | # TFs |
|----------:|-------------:|
| 1 |  328 |
| 2 |  161 |
| >=3 | 337 |