# QutRNA - minimal workflow

A minimal workflow for the local alignment of tRNAs reads with [parasail](https://github.com/jeffdaily/parasail/) of tRNAs from Oxford Nanopore direct RNA sequencing reads, adapted from [QutRNA](https://github.com/dieterich-lab/QutRNA).


```console
cd <LOCAL-DIRECTORY>
git clone https://github.com/wvdtoorn/QutRNA
```

## Alignment

The alignment workflow employs an optimal local alignment strategy using the implementation as provided by the [parasail](https://github.com/jeffdaily/parasail/) software.
In summary, optimal alignments may produce drastically different tRNA read mappings and are superior to heuristic alignments.

The strategy to assess the statistical significance is rooted in a simulation-based approach, which produces random alignments.
Briefly, we reverse input sequences.
Then, we compute alignments in **forward** orientation with all original reads and **reverse** orientation. We classify alignments in the **forward** orientation as true and the ones with **reverse** read orientation as false.

Alignment precision is then defined by TP/(TP+FP) and alignment recall by TP/(TP+FN) according to some score threshold *t*, where TP: true positive, FP: false positive, and FN: false negative.

We calculate an optimal threshold for a given precision and filter mapped reads accordingly.

We employ the alignment strategy separately on reads that are designated by Guppy as *fastq_pass* or *fastq_fail* and merge the results subsequently.

## Installation & Requirements:

Currently, the pipeline requires an X86\_64 architecture.

* [R 4.3](https://www.r-project.org/)
* Java >= 11
* [samtools](https://www.htslib.org/)
* [parasail v2.6.2](https://github.com/jeffdaily/parasail/archive/refs/tags/v2.6.2.tar.gz)
* (check `conda.yaml` in the repository)

We provide a YAML file to create a [conda](https://docs.conda.io/en/latest/) environment with all necessary software with the exception of [parasail](https://github.com/jeffdaily/parasail/). Unfortunatelly, no package for parasail exists in conda.

Create a conda environment with:
```console
conda env create -n qutrna -f <LOCAL-REPOSITORY>/conda.yaml
conda activate qutrna
```

### parasail

Follow [parasail compiling](https://github.com/jeffdaily/parasail?tab=readme-ov-file#autotools-build) instructions. 

If you install and compile from sources and you use conda, it is imperative to compile parasail within the conda environment. 

1. Create and activate conda environment (see above)
2. Download parasail [v2.6.2.tar.gz](https://github.com/jeffdaily/parasail/archive/refs/tags/v2.6.2.tar.gz)
3. Compile and install. Replace `<DESTINATION>` with your desired path for parasail, e.g. `~/miniconda3/envs/qutrna`:

```console
export PATH="$CONDA_PREFIX/x86_64-conda-linux-gnu/bin:$CONDA_PREFIX/x86_64-conda_cos6-linux-gnu/bin:$PATH"
conda install -y conda-forge::autoconf conda-forge::automake conda-forge::libtool conda-forge::gcc conda-forge::gxx conda-forge::zlib-devel
tar -zvpf v2.6.2.tar.gz
cd v2.6.2
autoreconf -fi
configure --prefix=<DESTINATION>
make
make install
```

Make sure, [parasail](https://github.com/jeffdaily/parasail) is installed in `$PATH`. 


### Snakemake workflow

The minimal workflow is implemented with [snakemake](https://github.com/snakemake/snakemake) and encompasses:

* tRNA alignment with [parasail](https://github.com/jeffdaily/parasail), 

The workflow can be configured with a single YAML file:
* `minimal_workflow/minimal_config.yaml` : analysis specific config, e.g.: parameters of tools.


### Sample table

Sample description `sample_table.tsv` must be TAB-separated file and contain the following columns:


| condition | sample_name | fastq |
| --------- | ----------- | ---------- |
| ...       | ...         | ...        |

`fastq`:  GZIPPED FASTQ sequencing reads is expected.

### Executing the workflow

Setup `minimal_workflow/minimal_config.yaml`, and `sample_table.tsv`.
Replace `<QUTRNA>` with the directory where you cloned the repository, add paths to the YAML files and start the workflow with:

```console
cd <QUTRNA>/minimal_workflow
snakemake -c 1 --configfile=minimal_config.yaml
```


#### Output

The output of the pipeline can  be found in the `<output_dir>` that you provide in `minimal_config.yaml`.

##### Alignment

Filtered BAM alignment files: `<output_dir>/resutls/bams/filtered/<sample>/<subsample>/`:

* `<output_dir>/resutls/bams/filtered/<sample>/<subsample>/<base-calling>_cutoff.txt`: Alignment score cutoff.
* `<output_dir>/resutls/bams/filtered/<sample>/<subsample>/orient~(fwd|rev)/<base-calling>_score.txt`: All alignment scores.
* `<output_dir>/resutls/bams/filtered/<sample>/<subsample>/orient~(fwd|rev)/<base-calling>.sorted.bam`: Unfiltered alignments.

Final BAM alignment files: `<output_dir>/resutls/bams/final/<bam-file>`.

# How to cite

Sun Y, Piechotta M, Naarmann-de Vries I, Dieterich C, Ehrenhofer-Murray AE. Detection of queuosine and queuosine precursors in tRNAs by direct RNA sequencing. Nucleic Acids Res. 2023 Nov 10;51(20):11197-11212. doi: 10.1093/nar/gkad826. PMID: 37811872; PMCID: PMC10639084.

# License

See LICENSE for details

# References

* [RNA modification mapping with JACUSA2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02676-0)

## Software
* [parasail](https://github.com/jeffdaily/parasail/)
