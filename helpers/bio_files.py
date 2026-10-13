def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Extracts the top Description per BLAST query and writes sorted results.

    For every 'Sequences producing significant alignments:' block the first
    data row's Description column is taken. Collected names are sorted
    alphabetically and written one per line.

    Arguments:
    input_file: str - path to the BLAST text report.
    output_file: str - path to the output txt file.

    Returns None.
    """
    descriptions: list[str] = []
    in_alignment_block = False

    with open(input_file, "r") as f:
        for line in f:
            stripped = line.strip()

            if stripped.startswith("Sequences producing significant alignments:"):
                in_alignment_block = True
                continue

            if in_alignment_block:
                if not stripped:
                    in_alignment_block = False
                    continue
                if stripped.lower().startswith("description"):
                    continue
                if stripped.startswith(("Score", "E value")):
                    continue
                parts = stripped.split("  ")
                description = parts[0].strip()
                if description:
                    descriptions.append(description)
                in_alignment_block = False

    descriptions = sorted(set(descriptions))

    with open(output_file, "w") as f:
        f.writelines(name + "\n" for name in descriptions)


import re


def parse_gbk_genes(input_gbk: str) -> list[dict]:
    """
    Parses a GenBank file and returns a list of CDS features.

    Each feature is a dictionary with keys 'name' and 'translation'.

    Arguments:
    input_gbk: str - path to the GenBank file.

    Returns a list of dicts.
    """
    with open(input_gbk, "r") as f:
        text = f.read()

    genes = []
    blocks = re.split(r"\n     CDS", text)
    for block in blocks[1:]:
        gene_match = re.search(r'/gene="([^"]+)"', block)
        trans_match = re.search(r'/translation="([^"]+)"', block)
        if gene_match and trans_match:
            translation = re.sub(r"\s+", "", trans_match.group(1))
            genes.append(
                {
                    "name": gene_match.group(1),
                    "translation": translation,
                }
            )
    return genes


def parse_gbk_genes(input_gbk: str) -> list[dict]:
    """
    Parses a GenBank file and returns a list of CDS features.

    Each feature is a dictionary with keys 'name' and 'translation'.

    Arguments:
    input_gbk: str - path to the GenBank file.

    Returns a list of dicts.
    """
    with open(input_gbk, "r") as f:
        text = f.read()

    genes = []
    blocks = re.split(r"\n     CDS", text)
    for block in blocks[1:]:
        gene_match = re.search(r'/gene="([^"]+)"', block)
        trans_match = re.search(r'/translation="([^"]+)"', block)
        if gene_match and trans_match:
            translation = re.sub(r"\s+", "", trans_match.group(1))
            genes.append(
                {
                    "name": gene_match.group(1),
                    "translation": translation,
                }
            )
    return genes


def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes,
    n_before: int = 1,
    n_after: int = 1,
    output_fasta: str = "neighbors.fasta",
) -> None:
    """
    Writes translations of genes flanking the genes of interest to a FASTA file.

    Arguments:
    input_gbk: str - path to the GenBank file.
    genes: str or list[str] - gene(s) of interest.
    n_before: int - number of genes before each gene of interest.
    n_after: int - number of genes after each gene of interest.
    output_fasta: str - path to the output FASTA file.

    Returns None.
    """
    if isinstance(genes, str):
        genes = [genes]
    genes_set = set(genes)

    all_genes = parse_gbk_genes(input_gbk)

    seen = set()
    selected = []

    for i, gene in enumerate(all_genes):
        if gene["name"] not in genes_set:
            continue
        start = max(0, i - n_before)
        end = min(len(all_genes), i + n_after + 1)
        for j in range(start, end):
            if j == i:
                continue
            neighbor = all_genes[j]
            if neighbor["name"] in seen:
                continue
            seen.add(neighbor["name"])
            selected.append(neighbor)

    with open(output_fasta, "w") as f:
        f.writelines(f">{g['name']}\n{g['translation']}\n" for g in selected)
