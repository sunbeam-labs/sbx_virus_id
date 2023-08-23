# sbx_virus_id



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