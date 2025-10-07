import csv
import os
from typing import Generator, TextIO


def parse_fasta(f: TextIO) -> Generator[tuple[str, str], None, None]:
    header_str = ""
    seq_str = ""
    for line in f.readlines():
        line = line.strip()
        if line.startswith(">"):
            if header_str:
                yield header_str, seq_str
            header_str = line
            seq_str = ""
        else:
            seq_str += line
    if header_str:
        yield header_str, seq_str


def write_fasta(record: list[str], f: TextIO) -> None:
    f.write(f"{record[0]}\n")
    f.write(f"{record[1]}\n")


idx = snakemake.input.idx  # type: ignore
fa = snakemake.input.fa  # type: ignore
output_fp = snakemake.output[0]  # type: ignore
log_fp = snakemake.log[0]  # type: ignore
contigs = {}

with open(idx) as f_idx:
    rd = csv.reader(f_idx, delimiter="\t", quotechar='"')
    for row in rd:
        if row[0] != "*":
            if int(row[2]) >= 1000:
                contigs[row[0]] = 1
            else:
                contigs[row[0]] = 0

with open(log_fp, "w") as f_log:
    f_log.write(f"Contigs: {contigs}")

with open(fa) as f_fa, open(output_fp, "w") as f_out:
    for header, seq in parse_fasta(f_fa):
        contig_name = header.split(" ")[0]
        if contigs[contig_name]:
            write_fasta([contig_name, seq], f_out)
