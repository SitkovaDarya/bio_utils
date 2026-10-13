"""
Helper functions for FASTQ filtering.

Provides file I/O and per-read validation used by filter_fastq.
"""

from typing import TextIO

Bounds = int | float | tuple[float, float]


def open_fastq(input_fastq: str) -> TextIO:
    """
    Opens a FASTQ file for reading.
    Arguments:
    input_fastq: str - path to the input FASTQ file.
    Returns the open file object.
    """
    return open(input_fastq, "r")


def read_one_record(f: TextIO) -> tuple[str, str, str, str] | None:
    """
    Reads a single FASTQ record from an open file.

    Arguments:
    f: TextIO - an open FASTQ file object.

    Returns a tuple (name, sequence, plus, quality) or None at end of file.
    """
    name = f.readline()
    if not name:
        return None
    sequence = f.readline().rstrip("\n")
    plus = f.readline().rstrip("\n")
    quality = f.readline().rstrip("\n")
    return name.rstrip("\n"), sequence, plus, quality


def normalize_bounds(bounds, default):
    """
    Converts a single number or a pair into an inclusive (low, high) tuple.

    Arguments:
    bounds: Bounds - a number (upper bound) or a (low, high) pair.
    default: tuple[float, float] - value used when bounds is None.

    Returns a tuple (low, high).
    """
    if bounds is None:
        return default
    if isinstance(bounds, (int, float)):
        return (0, float(bounds))
    low, high = bounds
    return (float(low), float(high))


def is_gc_in_bounds(sequence: str, gc_bounds: tuple[float, float]) -> bool:
    """
    Checks whether GC content of a sequence lies within bounds (in %).

    Arguments:
    sequence: str - nucleotide sequence.
    gc_bounds: tuple[float, float] - inclusive (low, high) GC range in percent.

    Returns bool.
    """
    if not sequence:
        return False
    gc_count = sum(sequence.upper().count(base) for base in ("G", "C"))
    gc_percent = 100.0 * gc_count / len(sequence)
    low, high = gc_bounds
    return low <= gc_percent <= high


def is_length_in_bounds(sequence: str, length_bounds: tuple[float, float]) -> bool:
    """
    Checks whether sequence length lies within bounds (inclusive).

    Arguments:
    sequence: str - nucleotide sequence.
    length_bounds: tuple[float, float] - inclusive (low, high) length range.

    Returns bool.
    """
    low, high = length_bounds
    return low <= len(sequence) <= high


def is_quality_ok(quality: str, threshold: float) -> bool:
    """
    Checks whether mean Phred33 quality of a read is >= threshold.

    Arguments:
    quality: str - quality string from a FASTQ record.
    threshold: float - minimum acceptable mean quality.

    Returns bool. Empty quality strings are rejected.
    """
    if not quality:
        return False
    total = sum(ord(ch) - 33 for ch in quality)
    return (total / len(quality)) >= threshold


def passes_all_filters(
    sequence: str,
    quality: str,
    gc_bounds: tuple[float, float],
    length_bounds: tuple[float, float],
    quality_threshold: float,
) -> bool:
    """
    Applies every filter to a single read.

    Arguments:
    sequence: str - nucleotide sequence.
    quality: str - quality string.
    gc_bounds: tuple[float, float] - inclusive GC range in percent.
    length_bounds: tuple[float, float] - inclusive length range.
    quality_threshold: float - minimum mean Phred33 quality.

    Returns True if the read passes all filters.
    """
    return (
        is_gc_in_bounds(sequence, gc_bounds)
        and is_length_in_bounds(sequence, length_bounds)
        and is_quality_ok(quality, quality_threshold)
    )
