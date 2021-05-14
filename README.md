## Genome Databases

All genome databases resides in the shared directory `/datalus/genomes/`.

Directory structure under each genome folder, using Arabidopsis thaliana as an example:

    (base) pzhou@node12:Athaliana $ pwd
    /datalus/genomes/Athaliana
    (base) pzhou@node12:Athaliana $ tree ./

```
├── 10.fasta  # genome fasta sequence
├── 10.fasta.fai
├── 15_intervals
│   ├── 01.chrom.bed  # chromosome size in BED format
│   ├── 01.chrom.sizes  # chromosome size file
│   ├── 11.gap.bed  # gap (N) coordinates
│   ├── forward.bed  # sequence name mappings between original (raw) fasta records and renamed (new) fasta records
│   ├── forward.chain  # UCSC chain file to convert raw (old) coordinates to new coordinates
│   ├── reverse.bed  # sequence name mappings between renamed (new) fasta records and original (raw) fasta records
│   └── reverse.chain # UCSC chain file to convert new coordinates to raw (old) coordinates
├── 21_dbs  # genome index databases
│   ├── bismark
│   │   ├── Bisulfite_Genome
│   │   │   ├── CT_conversion
│   │   │   │   ├── BS_CT.1.bt2
│   │   │   │   ├── BS_CT.2.bt2
│   │   │   │   ├── BS_CT.3.bt2
│   │   │   │   ├── BS_CT.4.bt2
│   │   │   │   ├── BS_CT.rev.1.bt2
│   │   │   │   ├── BS_CT.rev.2.bt2
│   │   │   │   └── genome_mfa.CT_conversion.fa
│   │   │   └── GA_conversion
│   │   │       ├── BS_GA.1.bt2
│   │   │       ├── BS_GA.2.bt2
│   │   │       ├── BS_GA.3.bt2
│   │   │       ├── BS_GA.4.bt2
│   │   │       ├── BS_GA.rev.1.bt2
│   │   │       ├── BS_GA.rev.2.bt2
│   │   │       └── genome_mfa.GA_conversion.fa
│   │   └── db.fasta
│   ├── blastn
│   │   ├── db.ndb
│   │   ├── db.nhr
│   │   ├── db.nin
│   │   ├── db.not
│   │   ├── db.nsq
│   │   ├── db.ntf
│   │   └── db.nto
│   ├── blastp
│   │   ├── db.pdb
│   │   ├── db.phr
│   │   ├── db.pin
│   │   ├── db.pot
│   │   ├── db.psq
│   │   ├── db.ptf
│   │   └── db.pto
│   ├── blat
│   │   ├── db.2bit
│   │   └── db.2bit.tile11.ooc
│   ├── bwa
│   │   ├── db.amb
│   │   ├── db.ann
│   │   ├── db.bwt
│   │   ├── db.fasta
│   │   ├── db.pac
│   │   └── db.sa
│   ├── gatk
│   │   ├── db.dict
│   │   ├── db.fasta
│   │   └── db.fasta.fai
│   ├── hisat2
│   │   ├── db.1.ht2
│   │   ├── db.2.ht2
│   │   ├── db.3.ht2
│   │   ├── db.4.ht2
│   │   ├── db.5.ht2
│   │   ├── db.6.ht2
│   │   ├── db.7.ht2
│   │   ├── db.8.ht2
│   │   ├── db.exon
│   │   └── db.ss
│   ├── lastn
│   │   ├── db.bck
│   │   ├── db.des
│   │   ├── db.prj
│   │   ├── db.sds
│   │   ├── db.ssp
│   │   ├── db.suf
│   │   └── db.tis
│   ├── lastp
│   │   ├── db.bck
│   │   ├── db.des
│   │   ├── db.prj
│   │   ├── db.sds
│   │   ├── db.ssp
│   │   ├── db.suf
│   │   └── db.tis
│   ├── salmon
│   │   └── db
│   │       ├── duplicate_clusters.tsv
│   │       ├── hash.bin
│   │       ├── header.json
│   │       ├── indexing.log
│   │       ├── quasi_index.log
│   │       ├── refInfo.json
│   │       ├── rsd.bin
│   │       ├── sa.bin
│   │       ├── txpInfo.bin
│   │       └── versionInfo.json
│   ├── snpeff
│   │   ├── Athaliana
│   │   │   ├── genes.gtf
│   │   │   ├── sequences.fa
│   │   │   └── snpEffectPredictor.bin
│   │   └── snpEff.config
│   └── star
│       ├── chrLength.txt
│       ├── chrNameLength.txt
│       ├── chrName.txt
│       ├── chrStart.txt
│       ├── exonGeTrInfo.tab
│       ├── exonInfo.tab
│       ├── geneInfo.tab
│       ├── Genome
│       ├── genomeParameters.txt
│       ├── Log.out
│       ├── SA
│       ├── SAindex
│       ├── sjdbInfo.txt
│       ├── sjdbList.fromGTF.out.tab
│       ├── sjdbList.out.tab
│       └── transcriptInfo.tab
├── 50_annotation
│   ├── 10.aa.fasta  # amino acid sequences for all transcripts
│   ├── 10.bed
│   ├── 10.desc.tsv
│   ├── 10.gff  # gff3 file with all transcript isoforms
│   ├── 10.gff.db
│   ├── 10.gtf
│   ├── 10.nt.fasta  # CDS sequences for all transcripts
│   ├── 10.sqlite
│   ├── 10.tsv
│   ├── 15.aa.fasta  # amino acid sequences for longest transcript isoforms
│   ├── 15.bed
│   ├── 15.desc.tsv
│   ├── 15.gff  # gff3 file for longest transcript isoforms
│   ├── 15.gff.db
│   ├── 15.gtf
│   ├── 15.nt.fasta  # CDS sequencces for longest transcript isoforms
│   └── 15.tsv
├── 55.rds
└── raw  # raw files downloaded
    ├── Arabidopsis_thaliana.TAIR10.50.gff3.gz
    ├── Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
    ├── raw.fasta.gz -> Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
    └── raw.gff.gz -> Arabidopsis_thaliana.TAIR10.50.gff3.gz
```
