"""
snakemake --rerun-triggers mtime -n all_ml_snpeff
snakemake --dag  --rerun-triggers mtime -R all_ml_snpeff | dot -Tsvg > ../results/img/control/dag_snpeff_db.svg

snakemake --jobs 60 \
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
  --jn job_ml.{name}.{jobid}.sh \
  -R all_ml_snpeff
"""

GTF_FILE = "../data/genomes/annotation/mirang.gtf.gz"
GFF_FILE = "../data/genomes/annotation/mirang.gff3.gz"

rule all_ml_snpeff:
    input: 
      vcf = expand( "../results/genotyping/annotated/{vcf_pre}_ann.vcf.gz", vcf_pre = "mirang_filtered" ),
      maked_load = expand( "../results/mutation_load/snp_eff/by_ind/masked/{sample}_masked.bed.gz", sample = SAMPLES ),
      expressed_load = expand( "../results/mutation_load/snp_eff/by_ind/expressed/{sample}_expressed.bed.gz", sample = SAMPLES ),
      load_in_roh = expand("../results/mutation_load/snp_eff/by_ind/{load_type}_in_roh/{sample}_{load_type}_in_roh.bed.gz", sample = SAMPLES, load_type = ["masked", "expressed", "fixed"] ),
      load_anc = expand( "../results/mutation_load/snp_eff/by_ind/{load_type}_anc/{sample}_{load_type}_anc.bed.gz", sample = SAMPLES, load_type = ["expressed", "fixed"] )

rule download_gtf:
    output:
      gtf = GTF_FILE
    shell:
      """
      wget https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Mirounga_angustirostris/annotation_releases/current/GCF_021288785.2-RS_2023_03/GCF_021288785.2_ASM2128878v3_genomic.gtf.gz
      mv GCF_021288785.2_ASM2128878v3_genomic.gtf.gz {output.gtf}
      """

rule download_gff:
    output:
      gff = GFF_FILE
    shell:
      """
      wget https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Mirounga_angustirostris/annotation_releases/current/GCF_021288785.2-RS_2023_03/GCF_021288785.2_ASM2128878v3_genomic.gff.gz
      mv GCF_021288785.2_ASM2128878v3_genomic.gff.gz {output.gff}
      """

rule create_snpeff_config:
    output:
      conf = "../results/mutation_load/snp_eff/snpEff.config"
    shell:
      """
      echo "# Mirounga angustirostris genome, version GCF_021288785.2" > {output.conf}
      echo "mirang.genome : mirang" >> {output.conf}
      """

rule extract_cds:
    input:
      fa = "../data/genomes/mirang.fa",
      gff = "../data/genomes/annotation/mirang.gff3.gz",
    output:
      cds = "../results/mutation_load/snp_eff/data/mirang/cds.fa.gz"
    params:
      cds_prefix = "../results/mutation_load/snp_eff/data/mirang/"
    conda: "gff3toolkit"
    shell:
      """
      gff3_to_fasta \
        -g {input.gff} \
        -f {input.fa} \
        -st cds \
        -d complete \
        -o {params.cds_prefix}/mirang
      
      mv {params.cds_prefix}/mirang_cds.fa {params.cds_prefix}/cds.fa
      gzip {params.cds_prefix}/cds.fa
      """

rule extract_prot:
    input:
      fa = "../data/genomes/mirang.fa",
      gff = "../data/genomes/annotation/mirang.gff3.gz"
    output:
      pep = "../results/mutation_load/snp_eff/data/mirang/protein.fa.gz"
    params:
      pep_prefix = "../results/mutation_load/snp_eff/data/mirang"
    conda: "gff3toolkit"
    shell:
      """
      gff3_to_fasta \
        -g {input.gff} \
        -f {input.fa} \
        -st pep \
        -d complete \
        -o {params.pep_prefix}/mirang
      
      mv {params.pep_prefix}/mirang_pep.fa {params.pep_prefix}/protein.fa
      gzip {params.pep_prefix}/protein.fa
      """

rule create_snpeff_db:
    input:
      fa = "../data/genomes/mirang.fa",
      gtf = GTF_FILE,
      cds = "../results/mutation_load/snp_eff/data/mirang/cds.fa.gz",
      prot = "../results/mutation_load/snp_eff/data/mirang/protein.fa.gz",
      conf = "../results/mutation_load/snp_eff/snpEff.config"
    output:
      snp_fa = "../results/mutation_load/snp_eff/data/genomes/mirang.fa",
      snp_gff = "../results/mutation_load/snp_eff/data/mirang/genes.gtf.gz"
    params:
      snpeff_path = "../results/mutation_load/snp_eff"
    resources:
      mem_mb=25600
    container: c_ml
    shell:
      """
      mkdir -p {params.snpeff_path}/data/mirang {params.snpeff_path}/data/genomes
      cd {code_dir}/{params.snpeff_path}/data/mirang
      ln -s {code_dir}/{input.gtf} ./genes.gtf.gz
      cd {code_dir}/{params.snpeff_path}/data/genomes
      ln -s {code_dir}/{input.fa} ./mirang.fa
      cd {code_dir}/{params.snpeff_path}
      snpEff build -Xmx24G -c {code_dir}/{input.conf} -dataDir $(pwd)/data -gtf22 -v mirang
      """

rule snpeff_link_vcf:
    input:
      snpeff_gff = "../results/mutation_load/snp_eff/data/mirang/genes.gtf.gz",
      vcf = "../results/genotyping/filtered/{vcf_pre}.vcf.gz"
    output:
      vcf_ln = "../results/mutation_load/snp_eff/{vcf_pre}.vcf.gz"
    params:
      snpeff_path = "../results/mutation_load/snp_eff"
    shell:
      """
      cd {code_dir}/{params.snpeff_path}
      ln -s {code_dir}/{input.vcf} ./
      """

rule run_snpeff:
    input:
      snpeff_gff = "../results/mutation_load/snp_eff/data/mirang/genes.gtf.gz",
      vcf = "../results/mutation_load/snp_eff/{vcf_pre}.vcf.gz"
    output:
      snpef_vcf = "../results/genotyping/annotated/{vcf_pre}_ann.vcf",
      report = "../results/mutation_load/snp_eff/{vcf_pre}_stats.html"
    params:
      snpeff_path = "../results/mutation_load/snp_eff"
    resources:
      mem_mb=25600
    container: c_ml
    shell:
      """
      cd {code_dir}/{params.snpeff_path}
      snpEff ann -Xmx24G -stats {wildcards.vcf_pre}_stats.html \
          -no-downstream \
          -no-intergenic \
          -no-intron \
          -no-upstream \
          -no-utr \
          -v \
          mirang {wildcards.vcf_pre}.vcf.gz > {code_dir}/{output.snpef_vcf}
      """

rule bgzip_vcf:
    input:
      vcf = "../results/genotyping/annotated/{vcf_pre}_ann.vcf"
    output:
      vcf = "../results/genotyping/annotated/{vcf_pre}_ann.vcf.gz",
      vcf_idx = "../results/genotyping/annotated/{vcf_pre}_ann.vcf.gz.tbi"
    conda: "popgen_basics"
    shell:
      """
      bgzip {input.vcf}
      tabix -p vcf {output.vcf}
      """

# at this point, the alleles need to be swapped to ancestral alleles
# within ancestral_alleles.smk

rule filter_load:
    input:
      vcf = "../results/ancestral_allele/mirang_filtered_{spec}_ann_aa.vcf.gz" 
    output:
      vcf = "../results/mutation_load/snp_eff/load_subset/mirang_filtered_{spec}_load.vcf.gz"
    resources:
      mem_mb=25600
    container: c_ml
    shell:
      """
      zcat {input.vcf} | \
        SnpSift filter "((exists LOF[*].NUMTR ) | ( ANN[*].IMPACT='HIGH' ) ) " | \
        bgzip > {output.vcf}
      tabix -p vcf {output.vcf}
      """

rule sample_order_load:
    input:
      vcf = "../results/mutation_load/snp_eff/load_subset/mirang_filtered_{spec}_load.vcf.gz"
    output:
      sample_order = "../results/mutation_load/snp_eff/load_subset/mirang_filtered_{spec}_load_sample_order.pop"
    log: "logs/load_sample_order_{spec}.log"
    resources:
      mem_mb=25600
    container: c_ml
    shell:
      """
      ./sh/sample_order_from_vcf.sh {input.vcf} > {output.sample_order}
      """

rule filter_fixed_load:
    input:
      vcf = "../results/ancestral_allele/mirang_filtered_ann_aa.vcf.gz",
      idx = "../results/ancestral_allele/mirang_filtered_ann_aa.vcf.gz.tbi",
      bed = "../results/mutation_load/snp_eff/snp_tally/fixed_in_{spec}.bed.gz"
    output:
      vcf = "../results/mutation_load/snp_eff/load_subset/fixed/mirang_filtered_{spec}_fixed_load.vcf.gz"
    resources:
      mem_mb=25600
    container: c_ml
    shell:
      """
      bcftools view \
        -Ov {input.vcf} \
        --regions-file {input.bed} | \
        SnpSift filter "((exists LOF[*].NUMTR ) | ( ANN[*].IMPACT='HIGH' ) ) " | \
        bgzip > {output.vcf}
      
      tabix -p vcf {output.vcf}
      """

rule sample_order_fixed_load:
    input:
      vcf = "../results/mutation_load/snp_eff/load_subset/fixed/mirang_filtered_{spec}_fixed_load.vcf.gz",
    output:
      sample_order = "../results/mutation_load/snp_eff/load_subset/fixed/mirang_filtered_{spec}_fixed_load_sample_order.pop"
    log: "logs/load_fixed_sample_order_{spec}.log"
    resources:
      mem_mb=25600
    container: c_ml
    shell:
      """
      ./sh/sample_order_from_vcf.sh {input.vcf} > {output.sample_order}
      """

rule masked_load:
    input: 
      vcf = lambda wc: expand( "../results/mutation_load/snp_eff/load_subset/mirang_filtered_{spec}_load.vcf.gz", spec = get_spec_from_sample(wc.sample) ),
      sample_order = lambda wc: expand( "../results/mutation_load/snp_eff/load_subset/mirang_filtered_{spec}_load_sample_order.pop", spec = get_spec_from_sample(wc.sample) )
    output:
      bed = "../results/mutation_load/snp_eff/by_ind/masked/{sample}_masked.bed.gz"
    resources:
      mem_mb=25600
    container: c_ml
    shell:
      """
      SAMPLE_IDX=$(awk '{{if($1=="{wildcards.sample}"){{print NR - 1}} }}' {input.sample_order})
      # heterozygous (masked load)
      zcat {input.vcf} | \
        SnpSift filter "( isHet(GEN[${{SAMPLE_IDX}}].GT) )" | \
        grep -v "^##" | \
        awk -v OFS="\t" -v s="{wildcards.sample}" '{{if(NR==1){{ for (i=1; i<=NF; ++i) {{ if ($i ~ s) c=i }} }} {{print $1,$2,$2,$c}} }}' | \
        sed 's/POS\tPOS/FROM\tTO/' | \
        gzip > {output.bed}
      """

rule expressed_load:
    input: 
      vcf = lambda wc: expand( "../results/mutation_load/snp_eff/load_subset/mirang_filtered_{spec}_load.vcf.gz", spec = get_spec_from_sample(wc.sample) ),
      sample_order = lambda wc: expand( "../results/mutation_load/snp_eff/load_subset/mirang_filtered_{spec}_load_sample_order.pop", spec = get_spec_from_sample(wc.sample) )
    output:
      bed = "../results/mutation_load/snp_eff/by_ind/expressed/{sample}_expressed.bed.gz"
    resources:
      mem_mb=25600
    container: c_ml
    shell:
      """
      SAMPLE_IDX=$(awk '{{if($1=="{wildcards.sample}"){{print NR - 1}} }}' {input.sample_order})
      # homozygous for affected allele (expressed load)
      # REF is affected allele
      EXPR_LOAD_REF="( (ANN[*].ALLELE = REF) & (isRef(GEN[${{SAMPLE_IDX}}].GT)) )"
      # ALT is affected allele
      EXPR_LOAD_ALT="( ( ! (ANN[*].ALLELE = REF)) & ((isVariant(GEN[${{SAMPLE_IDX}}].GT) & (isHom(GEN[${{SAMPLE_IDX}}].GT))) ))"
      zcat {input.vcf} | \
        SnpSift filter "( ${{EXPR_LOAD_REF}} | ${{EXPR_LOAD_ALT}} )" | \
        grep -v "^##" | \
        awk -v OFS="\t" -v s="{wildcards.sample}" '{{if(NR==1){{ for (i=1; i<=NF; ++i) {{ if ($i ~ s) c=i }} }} {{print $1,$2,$2,$c}} }}' | \
        sed 's/POS\tPOS/FROM\tTO/' | \
        gzip > {output.bed}
      """

rule fixed_load:
    input: 
      vcf = lambda wc: expand( "../results/mutation_load/snp_eff/load_subset/fixed/mirang_filtered_{spec}_fixed_load.vcf.gz", spec = get_spec_from_sample(wc.sample) ),
      sample_order = lambda wc: expand( "../results/mutation_load/snp_eff/load_subset/fixed/mirang_filtered_{spec}_fixed_load_sample_order.pop", spec = get_spec_from_sample(wc.sample) )
    output:
      bed = "../results/mutation_load/snp_eff/by_ind/fixed/{sample}_fixed.bed.gz"
    resources:
      mem_mb=25600
    container: c_ml
    shell:
      """
      SAMPLE_IDX=$(awk '{{if($1=="{wildcards.sample}"){{print NR - 1}} }}' {input.sample_order})
      # homozygous for affected allele
      # REF is affected allele
      EXPR_LOAD_REF="( (ANN[*].ALLELE = REF) & (isRef(GEN[${{SAMPLE_IDX}}].GT)) )"
      # ALT is affected allele
      EXPR_LOAD_ALT="( ( ! (ANN[*].ALLELE = REF)) & ((isVariant(GEN[${{SAMPLE_IDX}}].GT) & (isHom(GEN[${{SAMPLE_IDX}}].GT))) ))"
      zcat {input.vcf} | \
        SnpSift filter "( ${{EXPR_LOAD_REF}} | ${{EXPR_LOAD_ALT}} )" | \
        grep -v "^##" | \
        awk -v OFS="\t" -v s="{wildcards.sample}" '{{if(NR==1){{ for (i=1; i<=NF; ++i) {{ if ($i ~ s) c=i }} }} {{print $1,$2,$2,$c}} }}' | \
        sed 's/POS\tPOS/FROM\tTO/' | \
        gzip > {output.bed}
      """

rule expressed_load_anc:
    input: 
      vcf = lambda wc: expand( "../results/mutation_load/snp_eff/load_subset/mirang_filtered_{spec}_load.vcf.gz", spec = get_spec_from_sample(wc.sample) ),
      sample_order = lambda wc: expand( "../results/mutation_load/snp_eff/load_subset/mirang_filtered_{spec}_load_sample_order.pop", spec = get_spec_from_sample(wc.sample) )
    output:
      bed = "../results/mutation_load/snp_eff/by_ind/expressed_anc/{sample}_expressed_anc.bed.gz"
    resources:
      mem_mb=25600
    container: c_ml
    shell:
      """
      SAMPLE_IDX=$(awk '{{if($1=="{wildcards.sample}"){{print NR - 1}} }}' {input.sample_order})
      # homozygous for derived allele (expressed load - ancestral def)
      zcat {input.vcf} | \
        SnpSift filter "((isVariant(GEN[${{SAMPLE_IDX}}].GT)) & (isHom(GEN[${{SAMPLE_IDX}}].GT)))" | \
        grep -v "^##" | \
        awk -v OFS="\t" -v s="{wildcards.sample}" '{{if(NR==1){{ for (i=1; i<=NF; ++i) {{ if ($i ~ s) c=i }} }} {{print $1,$2,$2,$c}} }}' | \
        sed 's/POS\tPOS/FROM\tTO/' | \
        gzip > {output.bed}
      """

rule fixed_load_anc:
    input: 
      vcf = lambda wc: expand( "../results/mutation_load/snp_eff/load_subset/fixed/mirang_filtered_{spec}_fixed_load.vcf.gz", spec = get_spec_from_sample(wc.sample) ),
      sample_order = lambda wc: expand( "../results/mutation_load/snp_eff/load_subset/fixed/mirang_filtered_{spec}_fixed_load_sample_order.pop", spec = get_spec_from_sample(wc.sample) )
    output:
      bed = "../results/mutation_load/snp_eff/by_ind/fixed_anc/{sample}_fixed_anc.bed.gz"
    resources:
      mem_mb=25600
    container: c_ml
    shell:
      """
      SAMPLE_IDX=$(awk '{{if($1=="{wildcards.sample}"){{print NR - 1}} }}' {input.sample_order})
      # homozygous for derived allele (expressed load - ancestral def)
      zcat {input.vcf} | \
        SnpSift filter "((isVariant(GEN[${{SAMPLE_IDX}}].GT)) & (isHom(GEN[${{SAMPLE_IDX}}].GT)))" | \
        grep -v "^##" | \
        awk -v OFS="\t" -v s="{wildcards.sample}" '{{if(NR==1){{ for (i=1; i<=NF; ++i) {{ if ($i ~ s) c=i }} }} {{print $1,$2,$2,$c}} }}' | \
        sed 's/POS\tPOS/FROM\tTO/' | \
        gzip > {output.bed}
      """

# logically, THIS SHOULD be NULL
rule masked_in_roh:
    input:
      masked_load = "../results/mutation_load/snp_eff/by_ind/masked/{sample}_masked.bed.gz",
      roh = "../results/roh/bcftools/snp_based/bed/max_callable/roh_max_{sample}_on_mirang.bed"
    output:
      masked_in_roh = "../results/mutation_load/snp_eff/by_ind/masked_in_roh/{sample}_masked_in_roh.bed.gz"
    conda: "popgen_basics"
    shell:
      """
      intersectBed -a {input.roh} -b {input.masked_load} > {output.masked_in_roh}
      """

rule expressed_in_roh:
    input:
      expressed_load = "../results/mutation_load/snp_eff/by_ind/expressed/{sample}_expressed.bed.gz",
      roh = "../results/roh/bcftools/snp_based/bed/max_callable/roh_max_{sample}_on_mirang.bed"
    output:
      expressed_in_roh = "../results/mutation_load/snp_eff/by_ind/expressed_in_roh/{sample}_expressed_in_roh.bed.gz"
    conda: "popgen_basics"
    shell:
      """
      intersectBed -a {input.roh} -b {input.expressed_load} > {output.expressed_in_roh}
      """

rule fixed_in_roh:
    input:
      fixed_load = "../results/mutation_load/snp_eff/by_ind/fixed/{sample}_fixed.bed.gz",
      roh = "../results/roh/bcftools/snp_based/bed/max_callable/roh_max_{sample}_on_mirang.bed"
    output:
      fixed_in_roh = "../results/mutation_load/snp_eff/by_ind/fixed_in_roh/{sample}_fixed_in_roh.bed.gz"
    conda: "popgen_basics"
    shell:
      """
      intersectBed -a {input.roh} -b {input.fixed_load} > {output.fixed_in_roh}
      """

rule expressed_anc_in_roh:
    input:
      expressed_load = "../results/mutation_load/snp_eff/by_ind/expressed_anc/{sample}_expressed_anc.bed.gz",
      roh = "../results/roh/bcftools/snp_based/bed/max_callable/roh_max_{sample}_on_mirang.bed"
    output:
      expressed_in_roh = "../results/mutation_load/snp_eff/by_ind/expressed_anc_in_roh/{sample}_expressed_anc_in_roh.bed.gz"
    conda: "popgen_basics"
    shell:
      """
      intersectBed -a {input.roh} -b {input.expressed_load} > {output.expressed_in_roh}
      """

rule fixed_anc_in_roh:
    input:
      fixed_load = "../results/mutation_load/snp_eff/by_ind/fixed_anc/{sample}_fixed_anc.bed.gz",
      roh = "../results/roh/bcftools/snp_based/bed/max_callable/roh_max_{sample}_on_mirang.bed"
    output:
      fixed_in_roh = "../results/mutation_load/snp_eff/by_ind/fixed_anc_in_roh/{sample}_fixed_anc_in_roh.bed.gz"
    conda: "popgen_basics"
    shell:
      """
      intersectBed -a {input.roh} -b {input.fixed_load} > {output.fixed_in_roh}
      """