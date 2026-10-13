"""
bio_utils - a small collection of bioinformatics utilities.
"""

import os  # noqa: I001

from helpers.fastq import (
    open_fastq,
    read_one_record,
    normalize_bounds,
    passes_all_filters,
    Bounds,
)

from helpers.bio_files import (
    parse_blast_output,
    select_genes_from_gbk_to_fasta,
)


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Bounds = (0, 100),
    length_bounds: Bounds = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """
    Filters FASTQ reads and writes survivors into the 'filtered' folder.

    Arguments:
    input_fastq: str - path to the input FASTQ file.
    output_fastq: str - name of the output file inside 'filtered/'.
    gc_bounds: Bounds - inclusive GC range in percent. Single number = upper bound.
    length_bounds: Bounds - inclusive length range. Single number = upper bound.
    quality_threshold: float - minimum mean Phred33 quality.

    Returns None.
    """
    gc_pair = normalize_bounds(gc_bounds, (0, 100))
    length_pair = normalize_bounds(length_bounds, (0, 2**32))

    os.makedirs("filtered", exist_ok=True)
    output_path = os.path.join("filtered", output_fastq)

    fin = open_fastq(input_fastq)
    with open(output_path, "w") as fout:
        while True:
            record = read_one_record(fin)
            if record is None:
                break
            name, sequence, plus, quality = record
            if passes_all_filters(
                sequence, quality, gc_pair, length_pair, quality_threshold
            ):
                fout.write(f"{name}\n{sequence}\n{plus}\n{quality}\n")
    fin.close()
