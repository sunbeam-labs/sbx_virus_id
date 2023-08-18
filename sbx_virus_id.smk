# -*- mode: Snakemake -*-
#
# Rules for running Cenote-Taker2 and other tools in the viral id pipeline

VIRUS_FP = Cfg["all"]["output_fp"] / "virus"
TARGET_VIRUS_ID = [
    VIRUS_FP / "blastx" / f"{sample}.btf" for sample in Samples.keys()
]


try:
    BENCHMARK_FP
except NameError:
    BENCHMARK_FP = Cfg["all"]["output_fp"] / "benchmarks"
try:
    LOG_FP
except NameError:
    LOG_FP = Cfg["all"]["output_fp"] / "logs"


def get_ext_path() -> Path:
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


rule all_virus_id:
    input:
        TARGET_VIRUS_ID,


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


rule install_cenote_taker2:
    output:
        VIRUS_FP / "cenote_taker2" / ".installed",
    benchmark:
        BENCHMARK_FP / "install_cenote_taker2.tsv"
    log:
        LOG_FP / "install_cenote_taker2.log",
    params:
        ext_fp=str(get_ext_path()),
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
        contigs=ASSEMBLY_FP / "virus_id_megahit" / "{sample}_asm" / "final.contigs.fa",
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
        run_script=str(get_ext_path() / "Cenote-Taker2" / "run_cenote-taker2.py"),
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
        cd {params.out_dir}
        mkdir -p {params.sample}
        cd {params.sample}
        python {params.run_script} -c {input.contigs} -r {params.sample} -m 32 -t 32 -p true -db virion --cenote-dbs {params.db_fp} 2>&1 | tee {log}
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


rule build_virus_index:
    input:
        VIRUS_FP / "cenote_taker2" / "{sample}.fasta",
    output:
        VIRUS_FP / "cenote_taker2" / "{sample}.fasta.1.bt2",
    conda:
        "envs/sbx_virus_id.yml"
    threads: Cfg["sbx_virus_id"]["bowtie2_build_threads"]
    shell:
        "bowtie2-build --threads {threads} -f {input} {input}"


rule align_virus_reads:
    input:
        r1=QC_FP / host_decontam_Q() / "{sample}_1.fastq.gz",
        r2=QC_FP / host_decontam_Q() / "{sample}_2.fastq.gz",
        index=VIRUS_FP / "cenote_taker2" / "{sample}.fasta.1.bt2",
    output:
        temp(VIRUS_FP / "alignments" / "{sample}.sam"),
    params:
        index=VIRUS_FP / "cenote_taker2" / "{sample}.fasta",
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
        target=VIRUS_FP / "cenote_taker2" / "{sample}.fasta",
    conda:
        "envs/sbx_virus_id.yml"
    shell:
        """
        samtools view -bT {params.target} {input} > {output.bam}
        samtools sort -o {output.sorted} {output.bam}
        samtools index {output.sorted} {output.bai}
        """


rule calculate_virus_coverage:
    input:
        bam=VIRUS_FP / "alignments" / "{sample}.sorted.bam",
        idx=VIRUS_FP / "alignments" / "{sample}.sorted.bam.bai",
    output:
        VIRUS_FP / "alignments" / "{sample}.genomecoverage.txt",
    conda:
        "envs/sbx_virus_id.yml"
    shell:
        """
        samtools view -b {input.bam} | genomeCoverageBed -ibam stdin | grep -v 'genome' | perl scripts/coverage_counter.pl > {output}
        """


rule filter_virus_coverage:
    input:
        VIRUS_FP / "cenote_taker2" / "{sample}.fasta",
        VIRUS_FP / "alignments" / "{sample}.genomecoverage.txt",
    output:
        VIRUS_FP / "final_contigs.fasta",
    shell:
        "scripts/filter_virus_coverage.py"


# Install blast db:
# conda create -n blast
# conda activate blast
# conda install -c bioconda blast
# mkdir refseq_select_prot/
# cd refseq_select_prot/
# perl `which update_blastdb.pl` --decompress refseq_select_prot


rule virus_blastx:
    """Run diamond blastx on untranslated genes against a target db and write to blast tabular format."""
    input:
        VIRUS_FP / "cenote_taker2" / "{sample}.fasta",
    output:
        VIRUS_FP / "blastx" / "{sample}.btf",
    benchmark:
        BENCHMARK_FP / "run_virus_blastx_{sample}.tsv"
    log:
        LOG_FP / "run_virus_blastx_{sample}.log",
    params:
        blastdb=Cfg["sbx_virus_id"]["blast_db"],
    threads: Cfg["sbx_virus_id"]["blastx_threads"]
    resources:
        mem_mb=24000,
        runtime=720,
    conda:
        "envs/sbx_virus_id.yml"
    shell:
        """
        if [ -s {input} ]; then
            blastx \
            -query {input} \
            -db {params.blastdb} \
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
