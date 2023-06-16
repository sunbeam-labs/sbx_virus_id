import csv
from sunbeamlib.parse import parse_fasta, write_fasta


with open(snakemake.input.summary) as f_summary, open(snakemake.input.contigs) as f_contigs, open(snakemake.output[0], "w") as f_out:
    dr = csv.DictReader(f_summary, delimiter="\t")
    cd = {}
    phages = ["phage", "siphoviridae", "conjugative transposon"]
    for line in dr:
        if (all([x not in line["ORGANISM_NAME"].lower() for x in phages]) or snakemake.params["include_phages"]) and int(line["NUM_HALLMARKS"]) > 0:
            cd[line["ORIGINAL_NAME"]] = 1

    for header_str, seq_str in parse_fasta(f_contigs):
        if header_str.split(" ")[0] in cd:
            write_fasta([header_str, seq_str], f_out)