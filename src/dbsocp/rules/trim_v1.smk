
rule trim:
    """Trim away 5' and possible 3' handles on read1 and trim possible 3' handles on 
    read2."""
    output:
        interleaved_fastq="trimmed.fastq"
    input:
        r1_fastq="reads.1.fastq.gz",
        r2_fastq="reads.2.fastq.gz",
    log: "trimmed.fastq.log"
    threads: workflow.cores - 1  # rule tag needs one thread
    params:
        five_prime = "XNNN" + config["h1"] + "N"*len(config["barcode"]) + config["h2"],
        trim_len = sum([len(config["h1"]), len(config["barcode"]), len(config["h2"])])
    shell:
        "cutadapt"
        " -g '{params.five_prime};min_overlap={params.trim_len}...{config[h3]};optional'"
        " -A {config[h3]}"
        " --pair-filter 'any'"
        " -e 0.2"
        " -j {threads}"
        " -m 25"
        " --interleaved"
        " -o {output.interleaved_fastq}"
        " {input.r1_fastq}"
        " {input.r2_fastq}"
        " > {log}"

rule extract_DBS:
    """Extract barcode sequence from read1 FASTQ"""
    output:
        fastq="barcodes.fastq.gz"
    input:
         fastq="reads.1.fastq.gz"
    log: "barcodes.fastq.gz.log"
    threads: 20
    params:
        extract_len = len(config["h1"]),
        dbs_len_max = len(config["barcode"]) + 1,
        dbs_len_min = len(config["barcode"]) - 1
    shell:
        "cutadapt"
        " -g 'XNNN{config[h1]};min_overlap={params.extract_len}...{config[h2]}'"
        " -e 0.2"
        " --discard-untrimmed"
        " -j {threads}"
        " -m {params.dbs_len_min}"
        " -M {params.dbs_len_max}"
        " --max-n 0"
        " -o {output.fastq}"
        " {input.fastq}"
        " > {log}"

