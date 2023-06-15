"""
Run the DBS-OCP pipeline

This is a small wrapper around Snakemake that sets some default parameters
"""
from importlib_resources import files, as_file
import logging
import subprocess
import sys
from typing import List

from snakemake.utils import available_cpu_count

logger = logging.getLogger(__name__)


class SnakemakeError(Exception):
    pass


def add_arguments(parser):
    parser.add_argument(
        "-c",
        "--cores",
        metavar="N",
        type=int,
        default=available_cpu_count(),
        help="Run on at most N CPU cores in parallel. "
        "Default: %(default)s (all available cores).",
    )

    # This argument will not capture any arguments due to nargs=-1. Instead
    # parse_known_args() is used in __main__.py to add any arguments not
    # captured here to snakemake_args.
    smk_args = parser.add_argument_group("snakemake arguments")
    smk_args.add_argument(
        "snakemake_args",
        nargs=-1,
        help="Arguments passed to snakemake. For info about snakemake options run "
        "'snakemake --help'.",
    )


def main(args):
    try:
        run(cores=args.cores, snakemake_args=args.snakemake_args)
    except subprocess.CalledProcessError as e:
        print(f"Error in snakemake invocation: {e}", file=sys.stderr)
        sys.exit(e.returncode)
    sys.exit(0)


def run(
    cores: int = 1,
    snakefile: str = "Snakefile",
    workdir=None,
    snakemake_args: List[str] = None,
):
    with as_file(files("dbsocp") / snakefile) as snakefile_path:
        cmd = ["snakemake", "-s", str(snakefile_path), "--cores", str(cores)]

        # Set defaults
        cmd += ["--printshellcmds"]

        if workdir is not None:
            cmd += ["--directory", str(workdir)]

        if snakemake_args is not None:
            cmd += snakemake_args

        logger.debug(f"Command: {' '.join(cmd)}")
        subprocess.check_call(cmd)
