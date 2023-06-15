from pathlib import Path
import pytest
import shutil
import subprocess
import sys

from dbsocp.cli.init import init
from dbsocp.cli.config import change_config
from dbsocp.cli.run import run

TESTDATA = Path("tests/testdata")
TESTDATA_READ1_V1 = TESTDATA / "reads.1.fastq.gz"
TESTDATA_READ2_V1 = TESTDATA / "reads.2.fastq.gz"
TESTDATA_READ1_V2 = TESTDATA / "reads_v2.1.fastq.gz"
TESTDATA_READ2_V2 = TESTDATA / "reads_v2.2.fastq.gz"
TESTDATA_INDEX_V2 = (TESTDATA / "index_v2.fastq.gz").absolute()
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
        "sinto --version",
        "preseq",
    ]
    for tool in tools:
        print(f"'$ {tool}'")
        subprocess.run(tool.split(" "), stderr=sys.stdout)


def test_init(tmp_path):
    init(tmp_path / "analysis", TESTDATA_READ1_V1)


def test_config(tmp_path):
    workdir = tmp_path / "analysis"
    init(workdir, TESTDATA_READ1_V1)
    change_config(workdir / CONFIG_NAME, [("reference", str(TESTDATA_REFERENCE))])


def test_dryrun(tmp_path):
    workdir = tmp_path / "analysis"
    init(workdir, TESTDATA_READ1_V1)
    change_config(workdir / CONFIG_NAME, [("reference", str(TESTDATA_REFERENCE))])
    run(cores=1, workdir=workdir, snakemake_args=["--dryrun"])


@pytest.fixture(scope="module")
def _workdir_v1(tmp_path_factory):
    """
    This runs the pipeline using default parameters
    """
    path = tmp_path_factory.mktemp(basename="analysis-") / "analysis"
    init(path, TESTDATA_READ1_V1)
    change_config(path / CONFIG_NAME, [("reference", str(TESTDATA_REFERENCE))])
    try:
        run(cores=1, workdir=path, snakemake_args=["--cores", "1", "--keep-going"])
    except subprocess.CalledProcessError:
        print_all_logs(path)
    return path


@pytest.fixture
def workdir_v1(_workdir_v1, tmp_path):
    """Make a fresh copy of the prepared analysis directory"""
    path = tmp_path / "analysis"
    shutil.copytree(_workdir_v1, path)
    return path


@pytest.fixture(scope="module")
def _workdir_v2(tmp_path_factory):
    """
    This runs the pipeline using default parameters
    """
    path = tmp_path_factory.mktemp(basename="analysis-") / "analysis"
    init(path, TESTDATA_READ1_V2)
    change_config(
        path / CONFIG_NAME,
        [
            ("reference", str(TESTDATA_REFERENCE)),
            ("index", str(TESTDATA_INDEX_V2)),
            ("h2", "CATGACCTCTTGGAACTGTC"),
        ],
    )
    try:
        run(cores=1, workdir=path, snakemake_args=["--cores", "1", "--keep-going"])
    except subprocess.CalledProcessError:
        print_all_logs(path)
    return path


@pytest.fixture
def workdir_v2(_workdir_v2, tmp_path):
    """Make a fresh copy of the prepared analysis directory"""
    path = tmp_path / "analysis"
    shutil.copytree(_workdir_v2, path)
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
def test_run_creates_v1(workdir_v1, file):
    """Check that all expected files are created"""
    if not (workdir_v1 / file).exists():
        print(f"File {file} not found! See logs below:")
        print_all_logs(workdir_v1)
        raise AssertionError


@pytest.mark.parametrize("file", expected_files)
def test_run_creates_v2(workdir_v2, file):
    """Check that all expected files are created"""
    if not (workdir_v2 / file).exists():
        print(f"File {file} not found! See logs below:")
        print_all_logs(workdir_v2)
        raise AssertionError
