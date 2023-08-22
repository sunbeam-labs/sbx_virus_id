import csv
from sunbeamlib.parse import parse_fasta, write_fasta

contigs = {}

with open(snakemake.input.idx) as f_idx:
    rd = csv.reader(f_idx, delimiter="\t", quotechar='"')
    for row in rd:
        if row[0] != "*":
            if int(row[2]) >= 1000:
                contigs[row[0]] = 1
            else:
                contigs[row[0]] = 0

with open(snakemake.log, "w") as f_log:
    f_log.write(f"Contigs: {contigs}")

with open(snakemake.input.fa) as f_fa, open(snakemake.output[0], "w") as f_out:
    for header, seq in parse_fasta(f_fa):
        contig_name = header.split(" ")[0]
        if contigs[contig_name]:
            write_fasta((contig_name, seq), f_out)