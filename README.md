# pataphix

pataphix is a command line tool to process and filter PhiX sequences in fastq files and obtain QC metrics from them using bowtie2 for mapping. Those sequences are rarely used and often thrown away at the very beginning of analysis workflows, but sometimes you can have them bundled with your sample sequences. This tool helps to filter them out, and computes some QC metrics from them.

You can choose to write new fastq files without the PhiX sequences, as well as keep the alignment results. 

The metrics computed will be based on the PhiX alignment to its reference. Any difference, whether mismatch or indel, will be considered as an error. The number and rate of errors will be computer globally, and for each read position, in order to identify if there are specific parts of the reads that are more erroneous.

# Basic usage

```bash
pataphix --R1 data/my_reads_R1.fastq.gz --R2 data/my_reads_R2.fastq.gz --outdir my_data/
```

This will give you 4 files in the directory indicated by `--outdir`:

- "Mapping_stderr.log": The bowtie2 log.
- "phix.tsv": a .tsv file that lists the number of reads (or read pair) going to the PhiX genome, and the percentage of the total that it represents.
- "phix_error_estimates.tsv":  A .tsv file with 5 columns: The read of origin (R1 or R2), the position in the read, the error rate, the number of errors, the coverage of the position.
- "phix_metrics.tsv": A .tsv file with the metric names in the first column, and the metric value in the second column.

Multithreading can be used for the mapping step using `--threads`.

If you have single end reads, you can just provide your single reads with the `--R1` options without using the `--R2` option.

If you wish to filter your fastq files and keep only the reads that do not match the PhiX genome, you can use `--filter`.

If you wish to keep the resulting alignments (the bam file), you can use `--keep-alignments`.

If you only want to estimate the proportion of PhiX reads without doing anything else, you can use `--only-estimates`, a much faster version which will only parse the output of bowtie2 to give the percentage and number of PhiX reads.

# Help

```bash

usage: pataphix [-h] --R1 R1 [--R2 R2] [--tmpdir TMPDIR] [--outdir OUTDIR] [--only-estimates] [--filter] [--keep-alignments] [--loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--threads THREADS]
                [--version]

Check PhiX sequences in FASTQ files.

options:
  -h, --help            show this help message and exit
  --R1 R1               Path to the R1 FASTQ file.
  --R2 R2               Path to the R2 FASTQ file, if applicable.
  --tmpdir TMPDIR       Temporary directory to store intermediate files like indexes or bam files.
  --outdir OUTDIR       Output directory for results.
  --only-estimates      Only estimate the PhiX content without running the full analysis.
  --filter              Filter the reads based on PhiX content. If set, reads that do not map to PhiX will be written to the output directory in a file named similarly as the input file(s).
  --keep-alignments     Keep the alignments on the PhiX genome in SAM format after mapping.
  --loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Set the logging level.
  --threads THREADS     Number of threads to use for the mapping process.
  --version             Show the version of pataphix.

```
