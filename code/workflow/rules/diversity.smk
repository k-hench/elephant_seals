"""
snakemake --rerun-triggers mtime -n all_diversity
"""

rule all_diversity:
    input: 
      het = expand( "../results/het/het_{spec}.tsv", spec = [ 'mirang', 'mirleo' ]),
      pi = expand( "../results/pi/mirang_pi_dxy_{part}.tsv.gz", part = GENOME_PARTITIONS ),
      snp_dens = expand( "../results/snp_density/snp_dens_{spec}.tsv.gz", spec = [ 'mirang', 'mirleo' ] )

rule he_by_ind:
    input:
      vcf = "../results/genotyping/autosome/mirang_filtered_mirang_autosome.vcf.gz",
      inds = "../results/pop/inds_{spec}.pop"
    output:
      het = "../results/het/het_{spec}.tsv"
    container: c_popgen
    resources:
      mem_mb=15360
    shell:
      """
      vcftools \
        --gzvcf {input.vcf} \
        --keep {input.inds} \
        --mac 1 \
        --het \
        --stdout > {output.het}
      """

rule convert_012_by_spec:
    input:
      vcf = "../results/genotyping/filtered/mirang_filtered_mirang.vcf.gz",
      inds = "../results/pop/inds_{spec}.pop"
    output:
      g012 = "../results/genotyping/012/mirang_filtered_{spec}_012.tsv.gz"
    container: c_popgen
    resources:
      mem_mb=15360
    shell:
      """
      vcftools \
        --gzvcf {input.vcf} \
        --keep {input.inds} \
        --mac 1 \
        --012 \
        --stdout | gzip > {output.g012}
      """

wildcard_constraints:
    part="[^_]*"

rule all_bp_geno:
    input:
      vcf = "../results/genotyping/filtered/partitions/mirang_all_bp_{part}_filtered.vcf.gz"
    output:
      geno = "../results/genotyping/geno/partitions/mirang_all_bp_{part}.geno.gz"
    benchmark:
      'benchmark/genotyping/vcf_2_geno_mirang_all_bp_{part}.tsv'
    resources:
      mem_mb=15360
    container: c_popgen
    shell:
      """
      parseVCF -i {input.vcf} | gzip > {output.geno}
      """

rule pop_labels_all:
    input:
      p_mirang = "../results/pop/inds_mirang.pop",
      p_mirleo = "../results/pop/inds_mirleo.pop"
    output:
      p_all = "../results/pop/inds_all_labeled.pop"
    shell:
      """
      awk '{{print $1"\tmirang"}}' {input.p_mirang}  > {output.p_all}
      awk '{{print $1"\tmirleo"}}' {input.p_mirleo} >> {output.p_all}
      """

rule pi_and_dxy:
    input:
      geno = "../results/genotyping/geno/partitions/mirang_all_bp_{part}.geno.gz",
      pops = "../results/pop/inds_all_labeled.pop"
    output:
      pi_wind = "../results/pi/mirang_pi_dxy_{part}.tsv.gz"
    benchmark:
      'benchmark/pi/windows_pi_dxy_{part}.tsv'
    resources:
      mem_mb=20480
    container: c_popgen
    threads: 3
    shell:
      """
      popgenWindows \
        -w 100000 -s 25000 \
        --popsFile {input.pops} \
        -p mirang \
        -p mirleo \
        -g {input.geno} \
        -o {output.pi_wind} \
        -f phased \
        --writeFailedWindows \
        -T 3
      """

rule snp_density:
    input:
      vcf = "../results/genotyping/filtered/mirang_filtered_{spec}.vcf.gz",
    output:
      dens = "../results/snp_density/snp_dens_{spec}.tsv.gz"
    container: c_popgen
    resources:
      mem_mb=15360
    shell:
      """
      vcftools \
        --gzvcf {input.vcf} \
        --SNPdensity 100000 \
        --stdout | gzip > {output.dens}
      """