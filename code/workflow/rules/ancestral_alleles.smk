"""
snakemake --rerun-triggers mtime -n all_anc_allele
snakemake --rerun-triggers mtime   --use-singularity --singularity-args "--bind $CDATA" --use-conda -c 1 all_anc_allele

snakemake --jobs 100 \
  --latency-wait 30 \
  -p \
  --default-resources mem_mb=51200 threads=1 \
  --use-singularity \
  --singularity-args "--bind $CDATA" \
  --use-conda \
  --rerun-triggers mtime \
  --cluster '
    qsub \
      -V -cwd \
      -P fair_share \
      -l idle=1 \
      -l si_flag=1 \
      -pe multislot {threads} \
      -l vf={resources.mem_mb}' \
  --jn job_aa.{name}.{jobid}.sh \
  -R all_anc_allele
"""

rule all_anc_allele:
    input: "../results/ancestral_allele/mirang_filtered_ann_aa.vcf.gz"

rule extract_ancestral_hals:
    input:
      hal_anc = "../results/cactus/scratch/carnivora_set/tmp/steps-output/Anc52.hal",
      hal_mir = "../results/cactus/scratch/carnivora_set/tmp/steps-output/Anc56.hal"
    output:
      hal_anc = "../results/ancestral_allele/Anc52.hal",
      hal_mir = "../results/ancestral_allele/Anc56.hal"
    shell:
      """
      cp {input.hal_anc} {output.hal_anc}
      cp {input.hal_mir} {output.hal_mir}
      """

rule merge_ancestral_hals:
    input:
      hal_anc = "../results/ancestral_allele/Anc52.hal",
      hal_mir = "../results/ancestral_allele/Anc56.hal"
    output:
      hal = "../results/ancestral_allele/mir_ancestral.hal",
      check = "../results/ancestral_allele/mir_ancestral_summary.txt"
    resources:
      mem_mb=15360
    container: c_cactus
    shell:
      """
      cp {input.hal_anc} {output.hal}
      halAppendSubtree {output.hal} {input.hal_mir} Anc56 Anc56 --merge --hdf5InMemory --hdf5InMemory
      halStats {output.hal} > {output.check}
      """

rule anc_allele_tsv:
    input:
      hal = "../results/ancestral_allele/mir_ancestral.hal"
    output:
      tsv = "../results/ancestral_allele/mirang_anc52_snps.tsv"
    resources:
      mem_mb=15360
    container: c_cactus
    shell:
      """
      halSnps {input.hal} mirang Anc52 --tsv {output.tsv}
      """

rule pack_anc_allele_tsv:
    input:
      tsv = "../results/ancestral_allele/mirang_anc52_snps.tsv"
    output:
      gz = "../results/ancestral_allele/mirang_anc52_snps.tsv.gz"
    resources:
      mem_mb=15360
    shell:
      """
      gzip {input.tsv}
      """

rule anc_tree:
    input:
      hal = "../results/ancestral_allele/mir_ancestral.hal"
    output:
      nwk = "../results/ancestral_allele/mirang_anc_tree.nwk"
    container: c_cactus
    shell:
      """
      halStats {input.hal} --tree  > {output.nwk}
      """

rule snp_pos_and_alleles:
    input:
      vcf = "../results/genotyping/filtered/mirang_filtered.vcf.gz"
    output:
      snps_tsv = "../results/ancestral_allele/snps_vcf.tsv.gz"
    shell:
      """
      zgrep -v "^##" {input.vcf} | cut -f 1,2,4,5 | gzip > {output.snps_tsv}
      """

rule determine_anc_ref:
    input:
      tsv = "../results/ancestral_allele/mirang_anc52_snps.tsv.gz",
      snps_tsv = "../results/ancestral_allele/snps_vcf.tsv.gz"
    output:
      anc_tsv = "../results/ancestral_allele/anc_allele_assignment.tsv.gz",
      anc_bed = "../results/ancestral_allele/anc_allele_assignment.bed",
      tex_miss = "../results/tab/ancestral_allele_mismatches.tex"
    log:
      "logs/r_ancestral_ref_proposal.log"
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/ancient_alleles_assignment.R 2> {log} 1> {log}
      """

rule pack_aa_bed:
    input:
      bed = "../results/ancestral_allele/anc_allele_assignment.bed"
    output:
      gzbed = "../results/ancestral_allele/anc_allele_assignment.bed.gz"
    conda: "popgen_basics"
    shell:
      """
      bgzip {input.bed}
      tabix -s 1 -b 2 -e 3 {output.gzbed}
      """

rule annotate_vcf:
  input:
    vcf = "../results/genotyping/annotated/mirang_filtered_ann.vcf.gz",
    bed = "../results/ancestral_allele/anc_allele_assignment.bed.gz"
  output:
    vcf = temp( "../results/ancestral_allele/mirang_filtered_ann_aa.vcf" )
  conda: "popgen_basics"
  shell:
    """
    zcat {input.vcf} | \
      vcf-annotate -a {input.bed} \
        -d key=INFO,ID=AA,Number=1,Type=String,Description='Ancestral Allele' \
        -c CHROM,FROM,TO,INFO/AA > {output.vcf}
    """

rule convert_vcf_alleles:
  input:
    vcf = "../results/ancestral_allele/mirang_filtered_ann_aa.vcf"
  output:
    gzvcf = "../results/ancestral_allele/mirang_filtered_ann_aa.vcf.gz"
  resources:
      mem_mb=25600
  container: c_jvar
  shell:
    """
    java -jar /opt/jvarkit/dist/jvarkit.jar \
      vcffilterjdk \
      -f script.js {input.vcf} | \
      bgzip > {output.gzvcf}
    """