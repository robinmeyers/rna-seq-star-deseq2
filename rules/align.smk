def get_fq(wildcards):
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna().to_dict()
    else:
        # yes trimming, use trimmed data
        if not is_single_end(**wildcards):
            # paired-end sample
            return {'fq1' : expand("trimmed/{sample}-{unit}.1.fastq.gz", **wildcards),
                    'fq2' : expand("trimmed/{sample}-{unit}.2.fastq.gz", **wildcards)}
        # single end sample
        return {"fq1" : "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)}
            

rule align:
    input: unpack(get_fq)
    output:
        # see STAR manual for additional output files
        "star/{sample}-{unit}/Aligned.out.bam",
        "star/{sample}-{unit}/ReadsPerGene.out.tab"
    log:
        "logs/star/{sample}-{unit}.log"
    params:
        # path to STAR reference genome index
        index=config["ref"]["index"],
        # optional parameters
        extra="--quantMode GeneCounts --sjdbGTFfile {} {}".format(
              config["ref"]["annotation"], config["params"]["star"])
    threads: config["threads"]["star"]
    wrapper:
        "0.49.0/bio/star/align"
