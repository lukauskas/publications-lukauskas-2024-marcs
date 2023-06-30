rule genomic_bins:
    input:
        chromsizes=DATA_DIR / 'raw' / 'genome' / 'hg38.filtered.chrom.sizes'
    output:
        bed_gz=temp(INTERIM_OUTPUT_DIR / 'genome' / 'windows.{bin_size}.bp')
    params:
        bin_size=lambda w: int(w.bin_size)
    wildcard_constraints:
        bin_size='[1-9][0-9]+'
    conda:
        '../../envs/bedstats.yaml'
    threads:
        1
    shell:
        """
        bedtools makewindows -g {input.chromsizes} -w {params.bin_size} -i srcwinnum | gzip > {output.bed_gz}
        """
    

rule genomic_bins_excluding_blacklist:
    input:
        regions=rules.genomic_bins.output.bed_gz,
        blacklist=DATA_DIR / 'raw' / 'blacklist' / 'hg38-blacklist.v2.bed.gz'
    output:
        bed_gz=OUTPUT_DIR / 'genome' / 'windows.{bin_size}bp.noblacklist.bed.gz'
    wildcard_constraints:
        bin_size='[1-9][0-9]+'
    conda:
        '../../envs/bedstats.yaml'
    shell:
        """
        bedtools intersect -a {input.regions} -b {input.blacklist} -v | gzip > {output.bed_gz}
        """