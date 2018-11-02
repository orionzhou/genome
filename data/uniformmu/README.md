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

  - Table 1. Summary of genes with at least one exonic UniformMu
    insertion sites.

<!-- end list -->

``` r
require(tidyverse)
require(knitr)
require(kableExtra)
options(knitr.kable.NA = '')
dirw = '~/projects/genomes/data/uniformmu' 
fi = file.path(dirw, '/15.uniformmu.exon.tsv')
ti = read_tsv(fi)
tt = ti %>% count(n_mu)
kable(tt, format = 'markdown', booktabs = T, 
    align = 'r',
    col.names = c('# exonic insertions', '# genes'),
    caption = "Table 1.  Alignment statistics.") #%>%
```

| \# exonic insertions | \# genes |
| -------------------: | -------: |
|                 \>=3 |     5699 |
|                    1 |     5784 |
|                    2 |     3347 |

``` r
  #column_spec(1, bold=T) %>%#, width_min='7cm')
```
