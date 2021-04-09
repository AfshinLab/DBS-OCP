from pathlib import Path
import pytest
import shutil
import subprocess
import sys

from dbsocp.cli.init import init
from dbsocp.cli.config import change_config
from dbsocp.cli.run import run, SnakemakeError

TESTDATA = Path("tests/testdata")
TESTDATA_READ1 = TESTDATA / "reads.1.fastq.gz"
TESTDATA_READ2 = TESTDATA / "reads.2.fastq.gz"
TESTDATA_REFERENCE = (TESTDATA / "ref.fasta").absolute()
CONFIG_NAME = "dbsocp.yaml"


def print_all_logs(workdir: Path):
    """Print out all snakemake logs in workdir if they exists"""
    for log in workdir.glob("**/*.log"):
        print("==> ", log.name, " <==")
        with log.open() as f:
            print(f.read())


def test_environment():
    tools = [
        "python --version",
        "dbsocp --version",
        "snakemake --version",
        "starcode --version",
        "cutadapt --version",
        "bwa",
        "samtools --version",
        "snaptools --version",
        "preseq",
    ]
    for tool in tools:
        print(f"'$ {tool}'")
        subprocess.run(tool.split(" "), stderr=sys.stdout)


def test_init(tmp_path):
    init(tmp_path / "analysis", TESTDATA_READ1)


def test_config(tmp_path):
    workdir = tmp_path / "analysis"
    init(workdir, TESTDATA_READ1)
    change_config(workdir / CONFIG_NAME, [("reference", str(TESTDATA_REFERENCE))])


@pytest.fixture(scope="module")
def _workdir(tmp_path_factory):
    """
    This runs the pipeline using default parameters
    """
    path = tmp_path_factory.mktemp(basename="analysis-") / "analysis"
    init(path, TESTDATA_READ1)
    change_config(path / CONFIG_NAME, [("reference", str(TESTDATA_REFERENCE))])
    try:
        run(workdir=path, keepgoing=True)
    except SnakemakeError:
        pass
    return path


@pytest.fixture
def workdir(_workdir, tmp_path):
    """Make a fresh copy of the prepared analysis directory"""
    path = tmp_path / "analysis"
    shutil.copytree(_workdir, path)
    return path


expected_files = [
    "barcodes.fastq.gz",
    "barcodes.clstr.gz",
    "trimmed.barcoded.1.fastq.gz",
    "trimmed.barcoded.2.fastq.gz",
    "trimmed.barcoded.1_fastqc.html",
    "trimmed.barcoded.2_fastqc.html",
    "mapped.bam",
    "mapped.bam.bai",
    "mapped.snap",
    "mapped.snap.qc",
    "mapped.nsort.bam",
    "fragments.tsv.gz",
    "fragments.tsv.gz.tbi",
    "fragments.merged.tsv.gz",
    "fragments.merged.tsv.gz.tbi",
    "preseq_c_curve.txt",
    "multiqc_report.html",
]


@pytest.mark.parametrize("file", expected_files)
def test_run_creates(workdir, file):
    """Check that all expected files are created"""
    if not (workdir / file).exists():
        print(f"File {file} not found! See logs below:")
        print_all_logs(workdir)
        raise AssertionError
