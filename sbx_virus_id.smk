# -*- mode: Snakemake -*-
#
# Rules for running Cenote-Taker2 and other tools in the viral id pipeline

VIRUS_FP = Cfg["all"]["output_fp"] / "virus"
TARGET_VIRUS_ID = [
    VIRUS_FP / "hisss" / "summary" / "all_plot_summary.txt",
    VIRUS_FP / "hisss" / "summary" / "all_align_summary.txt",
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


rule all_virus_id:
    input:
        TARGET_VIRUS_ID,


rule virus_id_megahit_paired:
    input:
        r1=QC_FP / "decontam" / "{sample}_1.fastq.gz",
        r2=QC_FP / "decontam" / "{sample}_2.fastq.gz",
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
        "megahit_env.yml"
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
        loc=str(get_ext_path()),
    resources:
        runtime=2400,
    shell:
        """
        cd {params.loc}
        if [ ! -d "Cenote-Taker2" ]; then
            git clone https://github.com/mtisza1/Cenote-Taker2.git
        else
            echo "Cenote-Taker2 directory already exists"
        fi

        if {{ conda env list | grep 'cenote-taker2_env'; }} >/dev/null 2>&1; then
            echo "Cenote-Taker2 env already exists"
        else
            conda create --name cenote-taker2_env --file cenote_taker2_env_explicit.txt
        fi

        CONDA_BASE=$(conda info --base)
        source $CONDA_BASE/etc/profile.d/conda.sh
        conda activate cenote-taker2_env

        pip install phanotate

        # with all the options (75GB). The PDB database (--hhPDB) takes about 2 hours to download.
        cd Cenote-Taker2/
        python update_ct2_databases.py --hmm True --protein True --rps True --taxdump True --hhCDD True --hhPFAM True --hhPDB True >> {log}

        touch {output}
        """


rule cenote_taker2:
    input:
        contigs=expand(
            ASSEMBLY_FP / "virus_id_megahit" / "{sample}_asm" / "final.contigs.fa",
            sample=Samples.keys(),
        ),
        install=VIRUS_FP / "cenote_taker2" / ".installed",
    output:
        VIRUS_FP / "cenote_taker2" / "final_combined_virus_sequences_{sample}.fasta",
    benchmark:
        BENCHMARK_FP / "cenote_taker2_{sample}.tsv"
    log:
        LOG_FP / "cenote_taker2_{sample}.log",
    params:
        run_script=str(get_ext_path() / "Cenote-Taker2" / "run_cenote-taker2.py"),
        out_dir=str(VIRUS_FP / "cenote_taker2"),
    shell:
        """
        CONDA_BASE=$(conda info --base)
        source $CONDA_BASE/etc/profile.d/conda.sh
        conda activate cenote-taker2_env
        python {params.run_script} -c {input.contigs} -r {params.out_dir} -m 32 -t 32 -p true -db virion 2>&1 | tee {log}
        """


# Install blast db:
# conda create -n blast
# conda activate blast
# conda install -c bioconda blast
# mkdir refseq_select_prot/
# cd refseq_select_prot/
# perl `which update_blastdb.pl` --decompress refseq_select_prot

# rule build_virus_diamond_db:
#    """Use diamond makedb to create any necessary db indeces that don't exist."""
#    input:
#        Cfg["sbx_virus_id"]["blast_db"],
#    output:
#        VIRUS_FP / "diamond" / ".prep",
#    benchmark:
#        BENCHMARK_FP / "build_virus_diamond_db.tsv"
#    log:
#        LOG_FP / "build_virus_diamond_db.log",
#    conda:
#        "sbx_virus_id.yml"
#    shell:
#        """
#        diamond prepdb -d {input} 2>&1 | tee {log}
#        """


rule virus_blastx:
    """Run diamond blastx on untranslated genes against a target db and write to blast tabular format."""
    input:
        genes=VIRUS_FP
        / "cenote_taker2"
        / "final_combined_virus_sequences_{sample}.fasta",
    output:
        VIRUS_FP / "blastx" / "{sample}.btf",
    benchmark:
        BENCHMARK_FP / "run_virus_blastx_{sample}.tsv"
    log:
        LOG_FP / "run_virus_blastx_{sample}.log",
    params:
        blastdb=Cfg["sbx_virus_id"]["blast_db"],
    threads: Cfg["sbx_virus_id"]["blastx_threads"]
    conda:
        "sbx_virus_id.yml"
    shell:
        """
        if [ -s {input.genes} ]; then
            blastx \
            -q {input.genes} \
            --db {params.blastdb} \
            --outfmt 6 \
            --threads {threads} \
            --evalue 1e-10 \
            --max-target-seqs 2475 \
            --out {output} \
            2>&1 | tee {log}
        else
            echo "Caught empty query" >> {log}
            touch {output}
        fi
        """


rule hisss:
    input:
        expand(VIRUS_FP / "blastx" / "{sample}.btf", sample=Samples.keys()),
    output:
        VIRUS_FP / "hisss" / "summary" / "all_plot_summary.txt",
        VIRUS_FP / "hisss" / "summary" / "all_align_summary.txt",
    benchmark:
        BENCHMARK_FP / "hisss.tsv"
    log:
        LOG_FP / "hisss.log",
    conda:
        "hisss_env.yml"
    shell:
        """
        echo "HEY"
        """
