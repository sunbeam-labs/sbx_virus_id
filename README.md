<img src="https://github.com/sunbeam-labs/sunbeam/blob/stable/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_virus_id



## Introduction

sbx_virus_id is a [sunbeam](https://github.com/sunbeam-labs/sunbeam) extension for identifying viruses in samples. This pipeline uses [MEGAHIT](https://github.com/voutcn/megahit) or [SPAdes](https://github.com/ablab/spades) for assembly of contigs and [Cenote-Taker2](https://github.com/mtisza1/Cenote-Taker2) or [Virsorter2](https://github.com/jiarong/VirSorter2) for viral identification.

### Installation

```
sunbeam extend https://github.com/sunbeam-labs/sbx_virus_id.git
```

# Installing blast dbs

Install blast db:

```
conda create -n blast
conda activate blast
conda install -c bioconda blast
mkdir refseq_select_prot/
cd refseq_select_prot/
perl `which update_blastdb.pl` --decompress refseq_select_prot
```

Install viral blast db:

```
conda stuff from above ^^^
mkdir viral_prot/ && cd viral_prot/
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz && gzip -d viral.1.protein.faa.gz
makeblastdb -in viral.1.protein.faa -parse_seqids -title "viral" -dbtype prot
```

## Running

Run with sunbeam on the target `all_virus_id`,

```
sunbeam run --profile /path/to/project/ all_virus_id
```

### Options for config.yml

  - blast_db: path to blast db (default: "")
  - blastx_threads: number of threads for running blastx (default: 4)
  - bowtie2_build_threads: number of threads for running bowtie2-build (default: 4)
  - cenote_taker2_db: path to cenote-taker2 db (default: "")
  - virsorter_db: path to virsorter2 db (default: "")
  - host_decontam: Whether to run host decontamination (default: False)
  - include_phages: Whether to include phages in the output (default: False)
  - use_spades: Whether to use SPAdes instead of MEGAHIT (default: False)
  - use_virsorter: Whether to use Virsorter2 instead of Cenote-Taker2 (default: False)

## Legacy Installation

```
git clone https://github.com/sunbeam-labs/sbx_virus_id.git extensions/sbx_virus_id
cd extensions/sbx_virus_id
cat config.yml >> /path/to/sunbeam_config.yml
```