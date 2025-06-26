# Define RES variable for this file
RES = config["qutrna"]["output_dir"]

rule samtools_index:
  input: "{prefix}.sorted.bam",
  output: "{prefix}.sorted.bam.bai",
  conda: "qutrna",
  resources:
    mem_mb=2000
  log: RES + "/logs/samtools/index/{prefix}.log",
  shell: """
    samtools index {input:q} 2> {log:q}
  """


rule samtools_stats:
  input: "{prefix}.bam",
  output: "{prefix}.stats",
  conda: "qutrna",
  resources:
    mem_mb=2000
  log: RES + "/logs/samtools/stats/{prefix}.log",
  shell: """
    samtools stats {input:q} > {output:q} 2> {log:q}
  """


rule samtools_coverage:
  input: "{prefix}.bam",
  output: "{prefix}_coverage.tsv",
  conda: "qutrna",
  resources:
    mem_mb=2000
  log: RES + "/logs/samtools/coverage/{prefix}.log",
  shell: """
    samtools coverage {input:q} > {output:q} 2> {log:q}
  """


rule samtools_read_count:
  input: "{prefix}.bam",
  output: "{prefix}_read_count.txt",
  conda: "qutrna",
  resources:
    mem_mb=2000
  log: RES + "/logs/samtools/read_count/{prefix}.log",
  shell: """
    samtools view -c {input:q} > {output:q} 2> {log:q}
  """


rule samtools_multimappers:
  input: "{prefix}.bam",
  output: "{prefix}_multimappers.tsv",
  conda: "qutrna",
  resources:
    mem_mb=2000
  log: RES + "/logs/samtools/multimappers/{prefix}.log",
  shell: """
    ( samtools view {input:q} | \
      cut -f1,3 | \
      sort -k1,1 | \
      python  {workflow.basedir}/count_multimapper.py ) > {output:q} 2> {log:q}
  """



# TODO
rule samtools_get_score:
  input: "{prefix}.sorted.bam"
  output: "{prefix}_score.tsv"
  conda: "qutrna",
  resources:
    mem_mb=2000
  log: RES + "/logs/samtools/get_score/{prefix}.log",
  shell: """
    (samtools view {input:q} | cut -f 5 | sort | uniq -c | awk '{{print $2 "\t" $1}}') > {output:q} 2> {log:q}
  """

# FIXME
rule samtools_filter_by_cutoff:
  input: bam="results/bams/filtered/sample~{SAMPLE}/orient~fwd/sample~{SAMPLE}.sorted.bam",
         cutoff="results/bams/filtered/sample~{SAMPLE}/sample~{SAMPLE}_cutoff.txt",
  output: "results/bams/final/sample~{SAMPLE}/sample~{SAMPLE}.sorted.bam",
  conda: "qutrna",
  resources:
    mem_mb=2000
  log: RES + "/logs/samtools/filter_by_cutoff/sample~{SAMPLE}/calls.log",
  shell: """
    samtools view -b -q `cat {input.cutoff:q}` \
      -o {output:q} {input.bam:q} 2> {log:q}
  """


def _samtools_merge_input(wildcards):
  fnames = []
  for row in TBL.loc[[wildcards.SAMPLE]].itertuples():
    if READS == "fastq":
        fname = f"results/bams/final/sample~{row.sample_name}/calls.sorted.bam"
    elif READS == "bam":
      if wildcards.bam_type == "final":
        filter_ = FILTERS_APPLIED[-1]
        fname = f"results/bams/preprocessed/{filter_}/sample~{row.sample_name}/calls.sorted.bam"
      else:
        fname = f"results/bams/preprocessed/{{bam_type}}/sample~{row.sample_name}/calls.sorted.bam"
    else:
      raise Exception()
    fnames.append(fname)

  return fnames


rule samtools_merge:
  input: _samtools_merge_input,
  output: "results/bams/{bam_type}/{SAMPLE}.sorted.bam",
  conda: "qutrna",
  resources:
    mem_mb=10000
  log: RES + "/logs/samtools/merge/{SAMPLE}/{bam_type}.log",
  shell: """
    samtools merge {output:q} {input:q} 2> {log:q}
  """


# TODO
if config["parasail"]["lines"] > 0:
  def _samtools_merge_reads_input(wildcards):
    split_reads = checkpoints.parasail_split_reads.get(**wildcards).output[0]
    fnames, = glob_wildcards(
        os.path.join(
          split_reads,
          "part_{fname}.fastq.gz"))

    output_dir = "results/bams/mapped/sample~{SAMPLE}/orient~{ORIENT}/calls_split"
    return expand(os.path.join(output_dir, "part_{fname}_raw.bam"),
                  fname=fnames,
                  allow_missing=True)


  rule samtools_merge_split_reads:
    input: bams=_samtools_merge_reads_input,
           fasta=REF_FASTA,
    output: bam="results/bams/mapped/sample~{SAMPLE}/orient~{ORIENT}/calls_raw.bam",
            bai="results/bams/mapped/sample~{SAMPLE}/orient~{ORIENT}/calls_raw.bam.bai",
    conda: "qutrna",
    resources:
      mem_mb=10000
    log: RES + "/logs/samtools/merge_split_reads/bams/mapped/sample~{SAMPLE}/orient~{ORIENT}/calls.sorted.bam",
    shell: """
      ( samtools merge - {input.bams} | \
          samtools sort -o {output.bam:q} /dev/stdin && samtools index {output.bam:q} ) 2> {log:q}
    """
