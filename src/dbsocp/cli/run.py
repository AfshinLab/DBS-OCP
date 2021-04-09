"""
Run the DBS-OPC pipeline

This is a small wrapper around Snakemake that sets some default parameters
"""
from importlib_resources import path as resource_path
import logging
import sys

from snakemake import snakemake
from snakemake.utils import available_cpu_count

logger = logging.getLogger(__name__)


class SnakemakeError(Exception):
    pass


def add_arguments(parser):
    parser.add_argument(
        '--dryrun', '-n', default=False, action='store_true',
        help='Do not execute anything'
    )
    parser.add_argument(
        '--cores', '--jobs', '-j', metavar='N', type=int, default=available_cpu_count(),
        help='Run on at most N CPU cores in parallel. Default: %(default)s (all '
             'available cores)'
    )
    parser.add_argument(
        '--keepgoing', '-k', default=False, action='store_true',
        help='If one job fails, finish the others.'
    )
    parser.add_argument(
        '--unlock', default=False, action='store_true',
        help='Remove a lock on the working directory.'
    )
    parser.add_argument(
        '--delete-all-output', default=False, action='store_true',
        help="Remove all files generated by the snakemake workflow. Use together with "
             "-n/--dry-run to list files without actually deleting anything. "
             "Write-protected files are not removed. Nevertheless, use with care! "
             "Default: %(default)s"
    )
    parser.add_argument(
        '--force-run', '-R', nargs='*',
        help="Force the re-execution or creation of the given rules or files. Use this "
             "option if you changed a rule and want to have all its output in your "
             "workflow updated. Default: %(default)s"
    )
    dags = parser.add_mutually_exclusive_group()
    dags.add_argument(
        "--dag", default=False, action='store_true',
        help="Print the dag in the graphviz dot language (requires graphviz to be "
             "installed). Default: %(default)s. To get output to pdf file, pipe output "
             "into dot as follows: dbsocp run --dag | dot -Tpdf > dag.pdf"
    )
    dags.add_argument(
        "--filegraph", default=False, action='store_true',
        help="Print the file graph showing input/output file from rules in the "
             "graphviz dot language (requires graphviz to be installed). Default: "
             "%(default)s. To get output to pdf file, pipe output into dot "
             "as follows: dbsocp run --filegraph | dot -Tpdf > filegraph.pdf"
    )
    parser.add_argument(
        "--snakemake-kws", nargs=2, metavar=("KEY", "VALUE"), action="append",
        default=[],
        help=f"Give additional snakemake arguments not yet added to {__name__}. See "
             f"the Snakemake API ("
             f"https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html#"
             f") for options."
    )
    parser.add_argument(
        'targets', nargs='*', default=[],
        help="File(s) to create. If omitted, the full pipeline is run."
    )


def main(args):
    try:
        run(dryrun=args.dryrun,
            cores=args.cores,
            keepgoing=args.keepgoing,
            unlock=args.unlock,
            delete_all_output=args.delete_all_output,
            force_run=args.force_run,
            printdag=args.dag,
            printfilegraph=args.filegraph,
            snake_kws=dict(args.snakemake_kws),
            targets=args.targets)
    except SnakemakeError:
        sys.exit(1)
    sys.exit(0)


def run(
    dryrun: bool = False,
    cores: int = 4,
    keepgoing: bool = False,
    unlock: bool = False,
    delete_all_output: bool = False,
    force_run=None,
    printdag: bool = False,
    printfilegraph: bool = False,
    targets=None,
    workdir=None,
    snakefile="Snakefile",
    snake_kws=None,
):
    snake_kws = {} if snake_kws is None else snake_kws
    # snakemake sets up its own logging, and this cannot be easily changed
    # (setting keep_logger=True crashes), so remove our own log handler
    # for now
    logger.root.handlers = []
    with resource_path('dbsocp', snakefile) as snakefile_path:
        success = snakemake(
            snakefile_path,
            snakemakepath='snakemake',
            dryrun=dryrun,
            cores=cores,
            keepgoing=keepgoing,
            unlock=unlock,
            delete_all_output=delete_all_output,
            forcerun=force_run,
            printshellcmds=True,
            printdag=printdag,
            printfilegraph=printfilegraph,
            targets=targets,
            workdir=workdir,
            use_conda=True,
            printreason=dryrun,
            **snake_kws
        )
    if not success:
        raise SnakemakeError()
