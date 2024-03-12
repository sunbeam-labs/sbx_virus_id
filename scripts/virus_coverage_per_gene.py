from typing import TextIO


# Step 1: Pileup Coverage Calculation
def parse_coverage(pileup_f: TextIO) -> dict[str, list[int]]:
    coverage_data = {}
    for line in pileup_f:
        contig, _, _, coverage, _ = line.strip().split("\t", 4)
        if contig not in coverage_data:
            coverage_data[contig] = []
        coverage_data[contig].append(int(coverage))

    return coverage_data


# Step 2: Collating with Blast Results
def collate_coverage_blast(
    coverage_data: dict[str, list[int]], btf: TextIO
) -> list[str]:
    final_report = []

    for line in btf.readlines()[1:]:
        contig_id, gene_id, start_index, end_index = line.strip().split(",")
        start_index, end_index = int(start_index), int(end_index)

        if contig_id in coverage_data:
            contig_coverage = coverage_data[contig_id]

            # Calculate coverage within the gene region
            gene_coverage = sum(contig_coverage[start_index - 1 : end_index]) / (
                end_index - start_index + 1
            )
            final_report.append(
                {
                    "Contig_ID": contig_id,
                    "Gene_ID": gene_id,
                    "Gene_Coverage": gene_coverage,
                }
            )

    return final_report


# Step 3: Writing to Output
def write_output(final_report: list[str], output_f: TextIO) -> None:
    output_f.write("Contig_ID\tGene_ID\tGene_Coverage\n")
    for entry in final_report:
        output_f.write(
            f"{entry['Contig_ID']}\t{entry['Gene_ID']}\t{entry['Gene_Coverage']}\n"
        )


coverage_data = parse_coverage(snakemake.input.mpileup)
final_report = collate_coverage_blast(coverage_data, snakemake.input.btf)
write_output(final_report, snakemake.output.tsv)
