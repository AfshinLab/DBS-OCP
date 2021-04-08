"""
Convert BAM to Fragment file (BED)

For info on Fragment files go to:
https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments
"""
from itertools import takewhile
import logging
from collections import Counter
import sys

import pysam
from xopen import xopen

from dbsocp.utils import tqdm, Summary

logger = logging.getLogger(__name__)


def add_arguments(parser):
    parser.add_argument(
        "input", help="Input BAM file"
    )
    parser.add_argument(
            "-o", "--output", help="Output Fragment file of BED type"
        )
    parser.add_argument(
        "-m", "--min-mapq", type=int, default=20,
        help="Minimum mapping quality. Default: %(default)s"
    )
    parser.add_argument(
        "--barcode-tag", default="CB",
        help="SAM tag used to store cell barcode. Default: %(default)s"
    )


def main(args):
    cache = Counter()
    summary = Summary()
    current_chrom = None
    with xopen(args.output, "w") as outfile:
        for read, mate in parse_pairs(args.input):
            summary["Pairs read"] += 1
            if skip_pair(read, mate, args.min_mapq):
                summary["Pairs skipped"] += 1
                continue

            # Get 5-prime coordinates adjusted for soft-clipping
            read_5prime = compute_five_prime_coords(read)
            mate_5prime = compute_five_prime_coords(mate)

            # Get Tn5 adjusted start and end
            start = min(read_5prime, mate_5prime) + 4
            end = max(read_5prime, mate_5prime) - 5

            barcode = read.get_tag(args.barcode_tag)
            chromosome = read.reference_name

            if chromosome != current_chrom:
                write_chromosome(outfile, cache, current_chrom, summary)
                current_chrom = chromosome
                cache.clear()

            cache[(start, end, barcode)] += 1

        write_chromosome(outfile, cache, current_chrom, summary)

    summary["Duplication rate"] = 100 * summary["Pairs duplicate"] / \
        (summary["Pairs read"] - summary["Pairs skipped"])

    summary.print_stats()


def skip_pair(read, mate, threshold):
    """
    Skip read pair if
    - Either is unmapped
    - Not proper pair
    - Either has low mapping quality
    - Mapped to different chromosomes
    - Both mapped to same strand
    """
    return (read.is_unmapped or mate.is_unmapped) or \
           (not read.is_proper_pair) or \
           (read.mapping_quality < threshold or read.mapping_quality < threshold) or \
           (read.reference_name != mate.reference_name) or \
           (read.is_reverse == mate.is_reverse)


def compute_five_prime_coords(read):
    """
    Computes the 5' position in chromosome coordinates
    Assumes that the read is mapped and the CIGAR exists.

    From CellRanger-ATAC:
    https://github.com/10XGenomics/cellranger-atac/blob/5753003c286de55c193f02e1d66c847cb635aa60/lib/python/tools/peaks.py#L71
    """
    cigar = read.cigartuples
    if read.is_reverse:
        # Add the suffix clip to the alignment end position
        suffix_clip = sum([x[1] for x in takewhile(lambda x: x[0] == 4, reversed(cigar))])  # noqa: E501
        return read.reference_end + suffix_clip
    else:
        # Subtract the prefix clip from the alignment position
        prefix_clip = sum([x[1] for x in takewhile(lambda x: x[0] == 4, cigar)])
        return read.reference_start - prefix_clip


def write_chromosome(file, cache, chromosome, summary):
    if chromosome is not None:
        lines = [(chromosome, *pos, count) for pos, count in cache.items()]
        for line in sorted(lines, key=lambda x: (x[1], x[2])):
            summary["Fragments written"] += 1
            summary["Pairs duplicate"] += line[-1] - 1
            print(*line, sep="\t", file=file)


def parse_pairs(bam_file: str):
    """
    Yield read pairs for all properly paired read pairs in the input file.
    """
    cache = dict()
    # Fix for https://github.com/pysam-developers/pysam/issues/939
    save = pysam.set_verbosity(0)
    with pysam.AlignmentFile(bam_file) as openin:

        if openin.header["HD"]["SO"] != "coordinate":
            sys.exit("ERROR: Input BAM must be coordinate sorted.")

        for read in tqdm(openin, desc="Reading BAM"):
            if read.query_name in cache:
                yield read, cache.pop(read.query_name)
            else:
                cache[read.query_name] = read
    cache.clear()
    pysam.set_verbosity(save)
