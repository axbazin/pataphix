# pataphix

pataphix is a command line tool to process and filter PhiX sequences in fastq files and obtain QC metrics from them. Those sequences are rarely used and often thrown away at the very beginning of analysis workflows, but sometimes you can have them bundled with your sample sequences. This tool helps to filter them out, and computes some QC metrics from them.

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


