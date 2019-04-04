Maize UniformMu mutants characterization
================

* Raw GFF3 file mapping to B73 AGPv3 were obtained from maizeGDB
  * 79,407 insertion sites and related stocks information
* Insertion sites were lifted over to AGPv4
  * 77,085 sites
* The overlap of these sites with 40,239 AGPv4 gene models (longest alternative transcript) were checked:
  * 32,889 sites overlap w. at least one genes
  * 42,738 sites do not overlap w. any genic region (i.e., intergenic)

[Table listing all UniformMu insertions overlapping w. an AGPv4 gene](15.mu.genic.tsv)

* 15,597 genes have at least one insertion sites (UTR, CDS, intron)
* 13,631 have at least one exon insertion(s)
  * 7,905 have at least two exon insertions.

The W22 syntenic gene and mapping relationship (one-to-one, one-to-many, etc.) between B and W were added.

Structural integrity of these B73 genes in the W22 genome were checked based on variants called by [whole genome comparison between B73 and W22](https://github.com/orionzhou/wgc/blob/master/Rmd/wgc.md).

The assignment to a transcription factor family based on [PlantTFDB](http://planttfdb.cbi.pku.edu.cn/index.php?sp=Zma) was added.

Expression in 10 tissues from B, P and W were added.

* [Final table](16.gene.mu.tsv) lists all 13,631 genes with at least one exon insertions, with columns:
  * `gid_B73`: AGPv4 gene ID
  * `n_mu`: number exonic UniformMu insertion sites
  * `mid`, `sids`: mutant ID and stock ID
  * `fam`, `fam_size`, `fam_idx`: whether this gene is TF, the assigned TF family and member index within the family
  * `gid_W22`, `map_type`: W22 ortholog ID and mapping relationship
  * `impact`, `eff`: impact and effect of the change from B73 to W22





