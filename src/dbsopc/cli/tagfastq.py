"""
Tag FASTQ/FASTA headers with corrected barcodes.

ABOUT:

First the raw barcodes FASTQ is parser to get a dictionary:

    raw_barcodes[<HEADER>] = <RAW_BARCODE>

Then the corrected barcodes CLSTR file from starcode are parsed to get a dictionary:

    corrected_barcodes[<RAW_BARCODE>] = <CORRECTED_BARCODE>

The for each read-pair in the input FASTQ(s) the corrected barcode is recovered and
 used to tag the read by including it in the header.

    <HEADER> ==> <RAW_BARCODE> ==> <CORRECTED_BARCODE>
"""

from contextlib import ExitStack
from itertools import islice
import logging
import sys

import dnaio

from dbsopc.utils import tqdm, Summary

logger = logging.getLogger(__name__)

IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "M": "AC",
    "K": "GT",
    "S": "CG",
    "W": "AT",
    "H": "ACT",
    "B": "CGT",
    "V": "ACG",
    "D": "AGT",
    "N": "ACGT"
}


def add_arguments(parser):
    parser.add_argument(
        "uncorrected_barcodes",
        help="FASTQ/FASTA for uncorrected barcodes."
    )
    parser.add_argument(
        "corrected_barcodes",
        help="File for error corrected barcodes. Each line is tab delimited with the "
             "entries: (1) Corrected barcode, (2) Total read count, (3) Comma "
             "delimited barcode sequences corrected to (1)."
    )
    parser.add_argument(
        "input1",
        help="Input FASTQ/FASTA file. Assumes to contain read1 if given with second "
             "input file. If only input1 is given, input is assumed to be an "
             "interleaved. If reading from stdin is requested use '-' as a "
             "placeholder."
    )
    parser.add_argument(
        "input2", nargs='?',
        help="Input FASTQ/FASTA for read2 for paired-end read. Leave empty if using "
             "interleaved."
    )
    parser.add_argument(
        "--output1", "--o1",
        help="Output FASTQ/FASTA file name for read1. If not specified the result is "
             "written to stdout as interleaved. If output1 given but not output2, "
             "output will be written as interleaved to output1."
    )
    parser.add_argument(
        "--output2", "--o2",
        help="Output FASTQ/FASTA name for read2. If not specified but --o1/--output1 "
             "given the result is written as interleaved."
    )
    parser.add_argument(
        "--min-count", default=2, type=int,
        help="Minimum read count for barcode to be included in output."
    )
    parser.add_argument(
        "-p", "--pattern-match",
        help="IUPAC barcode string to match against corrected barcodes. Non-matched "
             "barcodes will be removed."
    )


def main(args):
    run_tagfastq(
        uncorrected_barcodes=args.uncorrected_barcodes,
        corrected_barcodes=args.corrected_barcodes,
        input1=args.input1,
        input2=args.input2,
        output1=args.output1,
        output2=args.output2,
        min_count=args.min_count,
        pattern_match=args.pattern_match,
    )


def run_tagfastq(
        uncorrected_barcodes: str,
        corrected_barcodes: str,
        input1: str,
        input2: str,
        output1: str,
        output2: str,
        min_count: int,
        pattern_match: str,
):
    logger.info("Starting")
    summary = Summary()
    # Get the corrected barcodes and create a dictionary pointing each raw barcode to
    # its canonical sequence.
    template = [set(IUPAC[base]) for base in pattern_match] if pattern_match else []
    with open(corrected_barcodes, "r") as reader:
        corrected_barcodes = parse_corrected_barcodes(reader, summary, template,
                                                      min_count)

    in_interleaved = not input2
    logger.info(f"Input is {'interleaved' if in_interleaved else 'paired'} FASTQ.")

    # If no output1 is given output is sent to stdout
    if not output1:
        logger.info("Writing output to stdout.")
        output1 = sys.stdout.buffer
        output2 = None

    out_interleaved = not output2
    logger.info(f"Output is {'interleaved' if out_interleaved else 'paired'} FASTQ.")

    # Parse input FASTA/FASTQ for read1 and read2, uncorrected barcodes and write output
    with ExitStack() as stack:
        reader = stack.enter_context(
            dnaio.open(input1, file2=input2, interleaved=in_interleaved, mode="r")
        )
        writer = stack.enter_context(
            dnaio.open(output1, file2=output2, interleaved=out_interleaved, mode="w")
        )
        uncorrected_barcode_reader = stack.enter_context(
            BarcodeReader(uncorrected_barcodes)
        )

        for read1, read2 in tqdm(reader, desc="Read pairs processed", disable=False):
            # Header parsing
            summary["Read pairs read"] += 1

            uncorrected_barcode_seq = uncorrected_barcode_reader.get_barcode(read1.name)
            corrected_barcode_seq = corrected_barcodes.get(uncorrected_barcode_seq)

            # Check if barcode was found and update header with barcode info.
            if corrected_barcode_seq:
                name = read1.name.split(maxsplit=1)[0]
                read1.name = f"{corrected_barcode_seq}:{name}\tCB:Z:{corrected_barcode_seq}"
                read2.name = f"{corrected_barcode_seq}:{name}\tCB:Z:{corrected_barcode_seq}"
            else:
                summary["Reads missing barcode"] += 1
                continue

            # Write to out
            summary["Read pairs written"] += 1
            writer.write(read1, read2)

    summary.print_stats(__name__)

    logger.info("Finished")


def parse_corrected_barcodes(open_file, summary, template, min_count):
    """
    Parse starcode cluster output and return a dictionary with raw sequences pointing
    to a corrected canonical sequence
    :param open_file: starcode tabular output file.
    :param summary: Summary instance
    :param template: List of sets with allowed bases for each position.
    :param min_count: int. Include barcode with read count grater than min_count
    :return: dict: raw sequences pointing to a corrected canonical sequence.
    """
    corrected_barcodes = dict()
    for cluster in tqdm(open_file, desc="Clusters processed"):
        canonical_seq, size, cluster_seqs = cluster.strip().split("\t", maxsplit=3)
        summary["Corrected barcodes"] += 1
        summary["Reads with corrected barcodes"] += int(size)
        summary["Uncorrected barcodes"] += len(cluster_seqs.split(","))

        if int(size) <= min_count:
            summary["Barcodes skipped"] += 1
            continue

        if template and not match_template(canonical_seq, template):
            summary["Barcodes not matching pattern"] += 1
            summary["Reads with barcodes not matching pattern"] += int(size)
            continue

        corrected_barcodes.update(
            {raw_seq: canonical_seq for raw_seq in cluster_seqs.split(",")}
        )

    return corrected_barcodes


def match_template(sequence: str, template) -> bool:
    if len(sequence) != len(template):
        return False

    for base, accepted_bases in zip(sequence, template):
        if base not in accepted_bases:
            return False
    return True


class BarcodeReader:
    def __init__(self, filename):
        self._cache = dict()
        self._file = dnaio.open(filename, mode="r")
        self.barcodes = self.parse()

    def parse(self):
        for barcode in self._file:
            yield barcode.name, barcode.sequence

    def get_barcode(self, read_name, maxiter=10):
        if read_name in self._cache:
            return self._cache.pop(read_name)

        for barcode_read_name, barcode_sequence in islice(self.barcodes, maxiter):
            # If read_name in next pair then parser lines are synced --> drop cache.
            if read_name == barcode_read_name:
                self._cache.clear()
                return barcode_sequence

            self._cache[barcode_read_name] = barcode_sequence
        return None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        self._file.close()
