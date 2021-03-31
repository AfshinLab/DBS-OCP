"""
Convert BAM to Fragment file (BED)
See: https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments
"""
from tqdm import tqdm
import logging
from collections import Counter
import sys

import pysam
from xopen import xopen

logger = logging.getLogger(__name__)


def add_arguments(parser):
    parser.add_argument(
        "input", help="Input BAM file"
    )
    parser.add_argument(
            "-o", "--output", help="Output Fragment file of BED type"
        )
    parser.add_argument(
        "-m", "--min-mapq", type=int, default=20, help="Minimum mapping quality. Default: %(default)s"
    )
    parser.add_argument(
        "--barcode-tag", default="CB", help="SAM tag used to store cell barcode. Default: %(default)s"
    )


def main(args):
    cache = Counter()
    current_chrom = None
    with xopen(args.output, "w") as outfile:
        for read, mate in parse_pairs(args.input):
            if read.mapping_quality < args.min_mapq or mate.mapping_quality < args.min_mapq:
                continue

            # Get Tn5 adjusted positions
            start = min(mate.reference_start, read.reference_start) + 4
            end = max(mate.reference_end, read.reference_end) - 5

            barcode = read.get_tag(args.barcode_tag)
            chromosome = read.reference_name

            if chromosome != current_chrom:
                write_chromosome(outfile, cache, current_chrom)
                current_chrom = chromosome
                cache.clear()

            cache[(start, end, barcode)] += 1

        write_chromosome(outfile, cache, current_chrom)


def write_chromosome(file, cache, chromosome):
    if chromosome is not None:
        lines = [(chromosome, *pos, count) for pos, count in cache.items()]
        for line in sorted(lines, key=lambda x: (x[1], x[2])):
            print(*line, sep="\t", file=file)


def parse_pairs(bam_file: str):
    """
    Yield read pairs for all properly paired read pairs in the input file.
    """
    cache = dict()
    save = pysam.set_verbosity(0)  # Fix for https://github.com/pysam-developers/pysam/issues/939
    with pysam.AlignmentFile(bam_file) as openin:

        if openin.header["HD"]["SO"] != "coordinate":
            sys.exit("ERROR: Input BAM must be coordinate sorted.")

        for read in tqdm(openin, desc="Reading BAM"):
            if read.query_name in cache:
                yield read, cache.pop(read.query_name)
            else:
                if not read.is_unmapped and not read.mate_is_unmapped and read.is_proper_pair:
                    cache[read.query_name] = read
    cache.clear()
    pysam.set_verbosity(save)
