# Encode Datasets

This directory contains a cache of encode dataset metadata.

It has been downloaded using the following command:

```
wget 'https://www.encodeproject.org/metadata/?control_type%21=%2A&status=released&perturbed=false&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assembly=GRCh38&type=Experiment&files.processed=true&files.preferred_default=true' -O encode_metadata.$(date +%F).tsv
```

And gzipped the files (`gzip encode_metadata.$(date +%F).tsv`).

Pay attention to `files.preferred_default=true` at the end of the query. The rest of the pipeline depends on it!.

## Statistics of 2021-11-05 dataset

Number of lines:

```
> zcat encode_metadata.2021-11-05.tsv.gz | wc -l
333312
```

First 3 lines:

```
❯ zcat encode_metadata.2021-11-05.tsv.gz | head -n3
File accession  File format     File type       File format type        Output type     File assembly   Experiment accession    Assay   Donor(s)        Biosample term id       Biosample term name     Biosample type  Biosample organism      Biosample treatments    Biosample treatments amount     Biosample treatments duration   Biosample genetic modifications methods     Biosample genetic modifications categories      Biosample genetic modifications targets Biosample genetic modifications gene targets    Biosample genetic modifications site coordinates        Biosample genetic modifications zygosity        Experiment target       Library made from       Library depleted in     Library extraction method       Library lysis method        Library crosslinking method     Library strand specific Experiment date released        Project RBNS protein concentration      Library fragmentation method    Library size range      Biological replicate(s) Technical replicate(s)  Read length     Mapped read length      Run type        Paired end      Paired with     Index of        Derived fromSize    Lab     md5sum  dbxrefs File download URL       Genome annotation       Platform        Controlled by   File Status     s3_uri  File analysis title     File analysis status    Audit WARNING   Audit NOT_COMPLIANT     Audit ERROR
ENCFF489FBK     bed narrowPeak  bed     narrowPeak      peaks   GRCh38  ENCSR217SET     DNase-seq       /human-donors/ENCDO499CXJ/      EFO:0002083     SW480   cell line       Homo sapiens                                                                                    DNA                                             2016-05-09      ENCODE                     22_1                                                     /files/ENCFF180EJG/, /files/ENCFF117RRB/        1636124 ENCODE Processing Pipeline      7d5814fdb89fa056b706294d51a23cf3                https://www.encodeproject.org/files/ENCFF489FBK/@@download/ENCFF489FBK.bed.gz                           released        s3://encode-public/2020/11/06/3c5aa045-5f16-4527-9f84-f470b126779c/ENCFF489FBK.bed.gz       ENCODE4 v3.0.0-alpha.2 GRCh38   released        low spot score
ENCFF181YJP     bigBed narrowPeak       bigBed  narrowPeak      peaks   GRCh38  ENCSR217SET     DNase-seq       /human-donors/ENCDO499CXJ/      EFO:0002083     SW480   cell line       Homo sapiens                                                                                    DNA                                             2016-05-09      ENCODE             2
2_1                                                     /files/ENCFF489FBK/     2566259 ENCODE Processing Pipeline      94eb4ac6fe21e2e4eb5b8e2ae2d19cc1                https://www.encodeproject.org/files/ENCFF181YJP/@@download/ENCFF181YJP.bigBed                           released        s3://encode-public/2020/11/06/0114bd64-6c0a-46b3-b9d8-737cd904d203/ENCFF181YJP.bigBed       ENCODE4 v3.0.0-alpha.2 GRCh38   released        low spot score
```

Last 3 lines

```
❯ zcat encode_metadata.2021-11-05.tsv.gz | tail -n3
ENCFF573YZM     bed narrowPeak  bed     narrowPeak      pseudoreplicated peaks  GRCh38  ENCSR621EXD     Histone ChIP-seq        /human-donors/ENCDO686OPO/      CL:1001608      foreskin fibroblast     primary cell    Homo sapiens                                                                            H3K36me3-human  DNA                                        2013-07-31       Roadmap                         1       1_1                                                     /files/ENCFF544XTF/, /files/ENCFF926FKE/, /files/ENCFF356LFX/   4036614 ENCODE Processing Pipeline      5126838fba9bea21409d529f739b2196                https://www.encodeproject.org/files/ENCFF573YZM/@@download/ENCFF573YZM.bed.gz                      released s3://encode-public/2021/01/30/c679fb23-dbe0-4fcc-b0e0-d538473b35a9/ENCFF573YZM.bed.gz   ENCODE4 v1.7.0 GRCh38   released
ENCFF663ECD     bigWig  bigWig          signal p-value  GRCh38  ENCSR621EXD     Histone ChIP-seq        /human-donors/ENCDO686OPO/      CL:1001608      foreskin fibroblast     primary cell    Homo sapiens                                                                            H3K36me3-human  DNA                                             2013-07-31      Roadmap                             1       1_1                                                     /files/ENCFF544XTF/, /files/ENCFF926FKE/        1044711538      ENCODE Processing Pipeline      7df3dcefd44d7817a4b7eb86aee95ec1                https://www.encodeproject.org/files/ENCFF663ECD/@@download/ENCFF663ECD.bigWig                           released        s3://encode-public/2021/01/30/6f18d26a-84dc-4ef7-b8be-c0262a0828f8/ENCFF663ECD.bigWig       ENCODE4 v1.7.0 GRCh38   released
ENCFF735PFE     bigBed narrowPeak       bigBed  narrowPeak      pseudoreplicated peaks  GRCh38  ENCSR621EXD     Histone ChIP-seq        /human-donors/ENCDO686OPO/      CL:1001608      foreskin fibroblast     primary cell    Homo sapiens                                                                            H3K36me3-human  DNA                                2013-07-31       Roadmap                         1       1_1                                                     /files/ENCFF573YZM/     7235269 ENCODE Processing Pipeline      6bfb618ed4e772509027f2e42616a7ca                https://www.encodeproject.org/files/ENCFF735PFE/@@download/ENCFF735PFE.bigBed                           released        s3://encode-public/2021/01/30/63409155-e89b-43a2-95fe-444a8956e2dd/ENCFF735PFE.bigBed       ENCODE4 v1.7.0 GRCh38   released
```