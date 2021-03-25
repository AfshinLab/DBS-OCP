from pathlib import Path
import pytest
import shutil
import subprocess

from dbsopc.cli.init import init
from dbsopc.cli.config import change_config
from dbsopc.cli.run import run

TESTDATA = Path("tests/testdata")
TESTDATA_READ1 = TESTDATA / "reads.1.fastq.gz"
TESTDATA_READ2 = TESTDATA / "reads.2.fastq.gz"
TESTDATA_REFERENCE = (TESTDATA / "ref.fasta").absolute()
CONFIG_NAME = "dbsopc.yaml"


def test_environment():
    subprocess.run(["python",  "--version"])
    subprocess.run(["snakemake", "--version"])
    subprocess.run(["starcode", "--version"])
    subprocess.run(["cutadapt", "--version"])
    subprocess.run(["bwa"])
    subprocess.run(["samtools", "--version"])
    subprocess.run(["snaptools", "--version"])


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

    run(workdir=path)
    return path


@pytest.fixture
def workdir(_workdir, tmp_path):
    """Make a fresh copy of the prepared analysis directory"""
    path = tmp_path / "analysis"
    shutil.copytree(_workdir, path)
    return path


def test_run(workdir):
    """Check that all expected files are created"""
    expected_files = [
        "barcodes.fastq.gz",
        "barcodes.clstr.gz",
        "trimmed.barcoded.1.fastq.gz",
        "trimmed.barcoded.2.fastq.gz",
        "trimmed.barcoded.1_fastqc.html",
        "trimmed.barcoded.2_fastqc.html",
        "mapped.bam",
        "mapped.snap",
        "mapped.snap.qc",
    ]
    for file in expected_files:
        assert (workdir / file).exists()
