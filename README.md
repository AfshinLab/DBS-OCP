[![CI](https://github.com/AfshinLab/DBS-OCP/workflows/CI/badge.svg?branch=main)](https://github.com/AfshinLab/DBS-OCP/actions?query=branch%3Amain)

# DBS-OCP
DBS for Open Chromatin Profiling in single cells/nuclei, a.k.a scATAC-seq/snATAC-seq.

- [Install](#install)
- [Running analysis](#running-analysis)

## Install

1. Clone the DBS-OPC github
    ```
    git clone https://github.com/AfshinLab/DBS-OCP 
    ```

2. Install the base environment for `environment.yaml` file using conda:
    ```
    conda env create -n dbsocp -f environment.yaml 
    ```

3. Install the DBS-OPC package:
    ```
    pip install . 
    ```
 
To install on the Uppmax Bianca cluster follow these [instructions](doc/bianca_install.rst).
     

## Running analysis

Required inputs:
- Paired FASTQ files
- BWA indexed reference 

1. Setup analysis directory called `my_workdir`:
    ```
    dbsocp init --read1 /path/to/read.1.fastq.gz my_workdir 
    ```
   This will create and populate the directory with symlinks to the paired
    FASTQ. 
2. Update configs, here the path to the BWA indexed genome reference is added
. To view current configs jsut run `dbsocp config` 
    ```
    dbsocp configs -s reference /path/to/reference
    ```
3. To run the full analysis type the following.  
    ```
    dbsocp run
    ```
   To only run the preprocessing steps (barcode correction, read trimming
   , barcode assignment and QC) include the argument `preprocess`. 


## Development

Information about ATAC-seq analysis can be found in the [here](doc/knowledge.md).

### Setup

To setup the development environment run:
```
conda env create -n dbsocp-dev -f environment.yaml
```

To install the package in editable mode with additional development dependancies run:
```
pip install -e .[dev]
```

### Code style

Code style is enforced using `black`. To run the formatter:
```
black src/ tests/
```

### Linting

Linting is done using `flake8`. To run the linter:
```
flake8 src/ tests/
flake8 --select=W292 --filename '*.yaml,*.yml'
```

### Testing

Testing is done using `pytest`. To run the tests:
```
pytest -v tests
```
