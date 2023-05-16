# -*- mode: Snakemake -*-
#
# Rules for running Cenote-Taker2 and other tools in the viral id pipeline

VIRUS_FP = Cfg["all"]["output_fp"] / "virus"
TARGET_VIRUS_ID = [VIRUS_FP / "hisss" / "all_virus.txt"]


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
        VIRUS_FP / "cenote_taker2" / ".installed"
    params:
        loc=str(get_ext_path())
    shell:
        """
        cd {params.loc}
        if [ ! -d "Cenote-Taker2" ]; then
            git clone https://github.com/mtisza1/Cenote-Taker2.git
        else
            echo "Cenote-Taker2 directory already exists"
        fi

        cd Cenote-Taker2
        conda env create --file cenote-taker2_env.yml
        conda activate cenote-taker2_env
        pip install phanotate
        conda install -c conda-forge -c bioconda hhsuite last=1282

        # with all the options (75GB). The PDB database (--hhPDB) takes about 2 hours to download.
        python update_ct2_databases.py --hmm True --protein True --rps True --taxdump True --hhCDD True --hhPFAM True --hhPDB True

        touch {output}
        """


rule cenote_taker2:
    input:
        contigs=expand(ASSEMBLY_FP / "virus_id_megahit" / "{sample}_asm" / "final.contigs.fa", sample=Samples.keys()),
        install=VIRUS_FP / "cenote_taker2" / ".installed",
    output:
        VIRUS_FP / "hisss" / "all_virus.txt"
    shell:
        "python /path/to/run_cenote-taker2.py -c /path/to/contigs -r output_directory -m 32 -t 32 -p true -db virion"