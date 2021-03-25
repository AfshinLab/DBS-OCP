# DBS-OCP
DBS for Open Chromatin Profiling in single nucleii 

## Install

1. Clone the DBS-OPC github
    ```
    git clone https://github.com/AfshinLab/DBS-OCP 
    ```

2. Install the base environment for `environment.yaml` file using conda:
    ```
    conda env create -n dbsopc -f environment.yaml 
    ```

3. Install the DBS-OPC package:
    ```
    pip install . 
    ```
    For development the package can be installed in editable mode using `pip
     install -e .`. 
     
To install on the Uppmax Bianca cluster follow these [instructions](doc/bianca_install.rst).
     

## Running analysis

Required inputs:
- Paired FASTQ files
- BWA indexed reference 

1. Setup analysis directory called `my_workdir`:
    ```
    dbsopc init --read1 /path/to/read.1.fastq.gz my_workdir 
    ```
   This will create and populate the directory with symlinks to the paired
    FASTQ. 
2. Update configs, here the path to the BWA indexed genome reference is added
. To view current configs jsut run `dbsopc config` 
    ```
    dbsopc configs -s reference /path/to/reference
    ```
3. To run the full analysis type the following.  
    ```
    dbsopc run
    ```
   To only run the preprocessing steps (barcode correction, read trimming
   , barcode assignment and QC) include the argument `preprocess`. 

