# -*- mode: Snakemake -*-
#
# Rules for running Cenote-Taker2 and other tools in the viral id pipeline

VIRUS_FP = Cfg["all"]["output_fp"] / "virus"


try:
    BENCHMARK_FP
except NameError:
    BENCHMARK_FP = Cfg["all"]["output_fp"] / "benchmarks"
try:
    LOG_FP
except NameError:
    LOG_FP = Cfg["all"]["output_fp"] / "logs"


def get_virus_ext_path() -> Path:
    ext_path = Path(sunbeam_dir) / "extensions" / "sbx_virus_id"
    if ext_path.exists():
        return ext_path
    raise Error(
        "Filepath for virus_id not found, are you sure it's installed under extensions/sbx_virus_id?"
    )


def host_decontam_Q() -> str:
    if Cfg["sbx_virus_id"]["host_decontam"]:
        return "decontam"
    else:
        return "cleaned"


def virus_sorter_input() -> Path:
    if Cfg["sbx_virus_id"]["use_spades"]:
        return ASSEMBLY_FP / "virus_id_spades" / "{sample}" / "scaffolds.fasta"
    else:
        return ASSEMBLY_FP / "virus_id_megahit" / "{sample}_asm" / "final.contigs.fa"


def virus_sorter_output() -> Path:
    if Cfg["sbx_virus_id"]["use_virsorter"]:
        return VIRUS_FP / "virsorter" / "{sample}.fasta"
    else:
        return VIRUS_FP / "cenote_taker2" / "{sample}.fasta"


rule all_virus_id:
    input:
        expand(VIRUS_FP / "blastx" / "{sample}.btf", sample=Samples.keys()),
        VIRUS_FP / "summary" / "all_align_summary.txt",


rule virus_id_megahit_paired:
    input:
        r1=QC_FP / host_decontam_Q() / "{sample}_1.fastq.gz",
        r2=QC_FP / host_decontam_Q() / "{sample}_2.fastq.gz",
    output:
        ASSEMBLY_FP / "virus_id_megahit" / "{sample}_asm" / "final.contigs.fa",
    benchmark:
        BENCHMARK_FP / "virus_id_megahit_paired_{sample}.tsv"
    log:
        LOG_FP / "virus_id_megahit_paired_{sample}.log",
    params:
        out_fp=str(ASSEMBLY_FP / "virus_id_megahit" / "{sample}_asm"),
    threads: 4
    conda:
        "envs/megahit_env.yml"
    resources:
        mem_mb=20000,
        runtime=720,
    shell:
        """
        ## turn off bash strict mode
        set +o pipefail

        ## sometimes the error is due to lack of memory
        exitcode=0
        if [ -d {params.out_fp} ]
        then
            echo "Clearing previous megahit directory..." > {log}
            rm -rf {params.out_fp}
        fi
        megahit -t {threads} -1 {input.r1} -2 {input.r2} -o {params.out_fp} --continue 2>&1 {log} || exitcode=$?

        if [ $exitcode -eq 255 ]
        then
            touch {output}
            echo "Empty contigs" 2>&1 | tee {log}
        elif [ $exitcode -gt 1 ]
        then
            echo "Check your memory" 2>&1 | tee {log}
        fi
        """


rule virus_id_spades_paired:
    input:
        r1=QC_FP / host_decontam_Q() / "{sample}_1.fastq.gz",
        r2=QC_FP / host_decontam_Q() / "{sample}_2.fastq.gz",
    output:
        ASSEMBLY_FP / "virus_id_spades" / "{sample}" / "scaffolds.fasta",
    benchmark:
        BENCHMARK_FP / "virus_id_spades_paired_{sample}.tsv"
    log:
        LOG_FP / "virus_id_spades_paired_{sample}.log",
    params:
        out_fp=str(ASSEMBLY_FP / "virus_id_spades" / "{sample}"),
    threads: 4
    conda:
        "envs/spades_env.yml"
    resources:
        mem_mb=20000,
        runtime=720,
    shell:
        """
        spades.py -1 {input.r1} -2 {input.r2} -t {threads} -o {params.out_fp} 2>&1 | tee {log}
        """


rule install_cenote_taker2:
    output:
        VIRUS_FP / "cenote_taker2" / ".installed",
    benchmark:
        BENCHMARK_FP / "install_cenote_taker2.tsv"
    log:
        LOG_FP / "install_cenote_taker2.log",
    params:
        ext_fp=str(get_virus_ext_path()),
        db_fp=Cfg["sbx_virus_id"]["cenote_taker2_db"],
    resources:
        runtime=2400,
    conda:
        "envs/cenote_taker2_env.yml"
    shell:
        """
        cd {params.ext_fp}
        bash scripts/install_cenote_taker2.sh {params.ext_fp} {output} {log} {params.db_fp}
        """


rule cenote_taker2:
    input:
        contigs=virus_sorter_input(),
        install=VIRUS_FP / "cenote_taker2" / ".installed",
    output:
        VIRUS_FP / "cenote_taker2" / "{sample}" / "final.contigs.fasta",
        VIRUS_FP
        / "cenote_taker2"
        / "{sample}"
        / "{sample}"
        / "{sample}_CONTIG_SUMMARY.tsv",
    benchmark:
        BENCHMARK_FP / "cenote_taker2_{sample}.tsv"
    log:
        LOG_FP / "cenote_taker2_{sample}.log",
    params:
        run_script=str(get_virus_ext_path() / "Cenote-Taker2" / "run_cenote-taker2.py"),
        out_dir=str(VIRUS_FP / "cenote_taker2"),
        sample="{sample}",
        db_fp=Cfg["sbx_virus_id"]["cenote_taker2_db"],
    resources:
        mem_mb=24000,
        runtime=720,
    conda:
        "envs/cenote_taker2_env.yml"
    shell:
        """
        if [[ {params.sample} == *"."* ]]; then
            echo "CenoteTaker2 doesn't allow sample names with '.' in them"
            exit 1
        fi

        cd {params.out_dir}
        mkdir -p {params.sample}
        cd {params.sample}
        
        python {params.run_script} -c {input.contigs} -r {params.sample} -m 32 -t 32 -p true -db virion --cenote-dbs {params.db_fp} 2>&1 | tee {log}
        """


rule install_virsorter:
    output:
        VIRUS_FP / "virsorter" / ".installed",
    benchmark:
        BENCHMARK_FP / "install_virsorter.tsv"
    log:
        LOG_FP / "install_virsorter.log",
    params:
        db_fp=Cfg["sbx_virus_id"]["virsorter_db"],
    resources:
        runtime=2400,
    threads: 4
    conda:
        "envs/virsorter_env.yml"
    shell:
        """
        # First check if directory exists and has files
        if [ -d {params.db_fp} ] && [ "$(ls -A {params.db_fp})" ]; then
            echo "VirSorter database already installed"
            touch {output}
            exit 0
        fi

        echo "Installing VirSorter database"
        virsorter setup -d {params.db_fp} -j 4
        touch {output}
        """


rule virsorter:
    input:
        contigs=virus_sorter_input(),
        install=VIRUS_FP / "virsorter" / ".installed",
    output:
        combined_viral=VIRUS_FP / "virsorter" / "{sample}" / "final-viral-combined.fa",
        scores=VIRUS_FP / "virsorter" / "{sample}" / "final-viral-score.tsv",
        boundaries=VIRUS_FP / "virsorter" / "{sample}" / "final-viral-boundary.tsv",
    benchmark:
        BENCHMARK_FP / "virsorter_{sample}.tsv"
    log:
        LOG_FP / "virsorter_{sample}.log",
    params:
        out_dir=str(VIRUS_FP / "virsorter" / "{sample}"),
        db_fp=Cfg["sbx_virus_id"]["virsorter_db"],
    resources:
        mem_mb=24000,
        runtime=720,
    threads: 4
    conda:
        "envs/virsorter_env.yml"
    shell:
        """
        virsorter run -w {params.out_dir} -i {input.contigs} --min-length 1000 -j {threads} --db-dir {params.db_fp} all
        """


rule filter_cenote_contigs:
    input:
        contigs=VIRUS_FP / "cenote_taker2" / "{sample}" / "final.contigs.fasta",
        summary=VIRUS_FP
        / "cenote_taker2"
        / "{sample}"
        / "{sample}"
        / "{sample}_CONTIG_SUMMARY.tsv",
    output:
        VIRUS_FP / "cenote_taker2" / "{sample}.fasta",
    params:
        include_phages=Cfg["sbx_virus_id"]["include_phages"],
    script:
        "scripts/filter_cenote_contigs.py"


rule filter_virsorter_contigs:
    input:
        contigs=VIRUS_FP / "virsorter" / "{sample}" / "final-viral-combined.fa",
    output:
        VIRUS_FP / "virsorter" / "{sample}.fasta",
    script:
        "scripts/filter_virsorter_contigs.py"


rule build_virus_index:
    input:
        virus_sorter_output(),
    output:
        f"{virus_sorter_output()}.1.bt2",
    conda:
        "envs/sbx_virus_id.yml"
    threads: Cfg["sbx_virus_id"]["bowtie2_build_threads"]
    shell:
        "bowtie2-build --threads {threads} -f {input} {input}"


rule align_virus_reads:
    input:
        r1=QC_FP / host_decontam_Q() / "{sample}_1.fastq.gz",
        r2=QC_FP / host_decontam_Q() / "{sample}_2.fastq.gz",
        index=f"{virus_sorter_output()}.1.bt2",
    output:
        temp(VIRUS_FP / "alignments" / "{sample}.sam"),
    params:
        index=str(virus_sorter_output()),
    threads: 6
    conda:
        "envs/sbx_virus_id.yml"
    shell:
        "bowtie2 -q --local -t --very-sensitive-local --threads {threads} --no-mixed --no-discordant -x {params.index} -1 {input.r1} -2 {input.r2} -S {output}"


rule process_virus_alignment:
    input:
        VIRUS_FP / "alignments" / "{sample}.sam",
    output:
        bam=temp(VIRUS_FP / "alignments" / "{sample}.bam"),
        sorted=temp(VIRUS_FP / "alignments" / "{sample}.sorted.bam"),
        bai=temp(VIRUS_FP / "alignments" / "{sample}.sorted.bam.bai"),
    params:
        target=str(virus_sorter_output()),
    conda:
        "envs/sbx_virus_id.yml"
    shell:
        """
        samtools view -bT {params.target} {input} > {output.bam}
        samtools sort -o {output.sorted} {output.bam}
        samtools index {output.sorted} {output.bai}
        """


rule calculate_mapping_stats:
    input:
        bam=VIRUS_FP / "alignments" / "{sample}.sorted.bam",
        idx=VIRUS_FP / "alignments" / "{sample}.sorted.bam.bai",
    output:
        VIRUS_FP / "alignments" / "{sample}.sorted.idxstats.tsv",
    params:
        ext_fp=str(get_virus_ext_path()),
    conda:
        "envs/sbx_virus_id.yml"
    shell:
        """
        samtools idxstats {input.bam} > {output}
        """


rule filter_virus_coverage:
    input:
        fa=virus_sorter_output(),
        idx=VIRUS_FP / "alignments" / "{sample}.sorted.idxstats.tsv",
    output:
        VIRUS_FP / "final_{sample}_contigs.fasta",
    log:
        LOG_FP / "filter_virus_coverage_{sample}.log",
    script:
        "scripts/filter_virus_coverage.py"


rule virus_blastx:
    """Run blastx on untranslated genes against a target db and write to blast tabular format."""
    input:
        VIRUS_FP / "final_{sample}_contigs.fasta",
    output:
        VIRUS_FP / "blastx" / "{sample}.btf",
    benchmark:
        BENCHMARK_FP / "run_virus_blastx_{sample}.tsv"
    log:
        LOG_FP / "run_virus_blastx_{sample}.log",
    params:
        blast_db=Cfg["sbx_virus_id"]["blast_db"],
    threads: Cfg["sbx_virus_id"]["blastx_threads"]
    resources:
        mem_mb=24000,
        runtime=720,
    conda:
        "envs/sbx_virus_id.yml"
    shell:
        """
        if [ -s {input} ]; then
            export BLASTDB=$(dirname {params.blast_db})
            blastx \
            -query {input} \
            -db $(basename {params.blast_db}) \
            -outfmt 6 \
            -num_threads {threads} \
            -evalue 0.05 \
            -max_target_seqs 100 \
            -out {output} \
            2>&1 | tee {log}
        else
            echo "Caught empty query" >> {log}
            touch {output}
        fi
        """


rule calculate_coverage:
    input:
        bam=VIRUS_FP / "alignments" / "{sample}.sorted.bam",
        idx=VIRUS_FP / "alignments" / "{sample}.sorted.bam.bai",
    output:
        VIRUS_FP / "alignments" / "{sample}.genomecoverage.txt",
    params:
        ext_fp=str(get_virus_ext_path()),
    conda:
        "envs/sbx_virus_id.yml"
    shell:
        """
        samtools view -b {input.bam} | genomeCoverageBed -ibam stdin | grep -v 'genome'| perl {params.ext_fp}/scripts/coverage_counter.pl > {output}
        """


rule combine_coverage_stats:
    input:
        cov=VIRUS_FP / "alignments" / "{sample}.genomecoverage.txt",
        stats=VIRUS_FP / "alignments" / "{sample}.sorted.idxstats.tsv",
    output:
        VIRUS_FP / "alignments" / "{sample}.align.summary.txt",
    benchmark:
        BENCHMARK_FP / "combine_coverage_stats_{sample}.tsv"
    log:
        LOG_FP / "combine_coverage_stats_{sample}.log",
    params:
        ext_fp=str(get_virus_ext_path()),
    conda:
        "envs/r_env.yml"
    shell:
        """
        Rscript {params.ext_fp}/scripts/combine_coverage_stats.R {input.cov} {input.stats} {output} 2>&1 | tee {log}
        """


rule all_summary:
    input:
        expand(
            VIRUS_FP / "alignments" / "{sample}.align.summary.txt",
            sample=Samples.keys(),
        ),
    output:
        VIRUS_FP / "summary" / "all_align_summary.txt",
    shell:
        """
        echo -e "Sample\tAlignTarget\tFractionCoverage\tTargetLength\tMappedReads" > {output}
        cat {input} >> {output}
        """
