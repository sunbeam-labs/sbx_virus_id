from sunbeamlib.parse import parse_fasta, write_fasta

# This does nothing but I'm leaving it in case we want to add custom filtering here later
with open(snakemake.input.contigs) as f_contigs, open(
    snakemake.output[0], "w"
) as f_out:
    for header_str, seq_str in parse_fasta(f_contigs):
        write_fasta((header_str, seq_str), f_out)
