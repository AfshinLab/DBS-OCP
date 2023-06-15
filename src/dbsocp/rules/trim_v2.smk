"""
Trim V2 reads

Layout

    5' - h1 - barcode - h2 - x - insert - h3 - 3'
    R1                         |--->    
    R2                              <---|
    I2     <---------------|

x: is a part of h2 that is used as sequencing priming site
"""
from dbsocp.utils import revcomp

rule trim:
    """Trim possible 3' handles on read1 and read2."""
    output:
        interleaved_fastq="trimmed.fastq"
    input:
        r1_fastq="reads.1.fastq.gz",
        r2_fastq="reads.2.fastq.gz",
    log: "trimmed.fastq.log"
    threads: workflow.cores - 1  # rule tag needs one thread
    params:
        h3 = config["h3"]
    shell:
        "cutadapt"
        " -a {params.h3}"
        " -A {params.h3}"
        " --pair-filter 'any'"
        " -e 0.1"
        " -j {threads}"
        " -m 25"
        " --interleaved"
        " -o {output.interleaved_fastq}"
        " {input.r1_fastq}"
        " {input.r2_fastq}"
        " > {log}"


rule symlink_index:
    output:
        fastq = "index.fastq.gz"
    params:
        index = config["index"]
    shell:
        "ln -s {params.index} {output.fastq}"


rule extract_DBS:
    """Extract barcode sequence from index 2 FASTQ"""
    output:
        fastq = "barcodes.fastq.gz"
    input:
        fastq = "index.fastq.gz"
    log: "barcodes.fastq.gz.log"
    threads: 20
    params:
        h2_comp = revcomp(config["h2"]),
        h1_comp = revcomp(config["h1"]),
        dbs_len_max = len(config["barcode"]) + 1,
        dbs_len_min = len(config["barcode"]) - 1
    shell:
        "cutadapt"
        " -g '^{params.h2_comp}...{params.h1_comp}'"
        " -e 0.1"
        " --discard-untrimmed"
        " -j {threads}"
        " -m {params.dbs_len_min}"
        " -M {params.dbs_len_max}"
        " --max-n 0"
        " -o {output.fastq}"
        " {input.fastq}"
        " > {log}"
