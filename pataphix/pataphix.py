#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
pataphix - A Python utility for checking PhiX sequences in FASTQ files.
"""
# default libraries
import argparse
import sys
import subprocess
import logging
from importlib.metadata import distribution
from pathlib import Path
from tempfile import TemporaryDirectory, gettempdir
from typing import List, IO
import os

# third-party libraries
import pysam


def is_bowtie2_installed():
    """
    Tries to run bowtie2 to see if it is installed and available in this scope
    """
    try:
        subprocess.run(["bowtie2", "--version"], capture_output=True)
        logging.info("bowtie2 is installed. Proceeding.")
    except FileNotFoundError as e:
        raise Exception(
            "bowtie2 is not installed or an error is raised when using it. Command run: 'bowtie2 --version'"
        ) from e


def index_fasta(tmpdir: TemporaryDirectory[str]) -> str:
    """
    Generate the bowtie2 index of the PhiX fasta file included in the code for mapping.
    """
    fasta = Path(os.path.dirname(os.path.realpath(__file__)) + "/data/phix.fasta")
    bowtie2_cmd: List[str] = [
        "bowtie2-build",
        str(fasta),
        tmpdir.name + "/PhiX",
    ]
    logging.debug(f"Launching: {' '.join(bowtie2_cmd)}")
    subprocess.run(bowtie2_cmd, capture_output=True)
    logging.debug(f"Done indexing PhiX fasta file: {fasta}")
    logging.debug(f"Index files are in: {tmpdir.name}")
    logging.debug(f"Index files are: {os.listdir(tmpdir.name)}")
    return tmpdir.name + "/PhiX"


def map_reads(
    R1: Path,
    R2: Path | None,
    index: str,
    threads: int,
    outdir: Path,
    filter: bool = False,
) -> subprocess.Popen[str]:
    """
    Map the reads from the FASTQ files to the PhiX index using bowtie2.
    """

    bowtie2_cmd: List[str] = [
        "bowtie2",
        "--no-unal",
    ]

    if filter:
        filter_option = "--un-conc-gz" if R2 else "--un-gz"
        filter_path = os.path.basename(R1).replace(".fastq", "_phix_filtered.fastq")
        if R2 is not None:
            filter_path = (
                os.path.basename(R1).split("_R1_")[0] + "_phix_filtered.fastq.gz"
            )
        bowtie2_cmd.extend([filter_option, str(outdir / filter_path)])
    bowtie2_cmd.extend(
        [
            "--end-to-end",
            "--threads",
            str(threads),
            "-x",
            index,
        ]
    )

    if R2 is not None:
        bowtie2_cmd.extend(["-1", str(R1), "-2", str(R2)])
    else:
        bowtie2_cmd.extend(["-U", str(R1)])

    logging.debug(f"Launching: {' '.join(bowtie2_cmd)}")
    mapping_process = subprocess.Popen(
        bowtie2_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=1,
        universal_newlines=True,
    )
    return mapping_process


def deal_with_mapping(process_output: IO[str], outdir: Path, keep: bool = False):
    """
    Process the output of the bowtie2 mapping.
    If keep is True, the SAM file will be read; otherwise, the process output will be read.
    """
    # process_output is the stdout of the bowtie2 process
    # then wait for the process to finish and save the SAM file
    bam_file = outdir / "mapped_phix.bam"
    instream = pysam.AlignmentFile(process_output, "r")
    if keep:
        bamfout = pysam.AlignmentFile(bam_file, "wb", template=instream)
    for line in instream:
        # do stuff with the reads
        if keep:
            bamfout.write(line)
    if keep:
        bamfout.close()
        logging.info(f"Mapping results saved to {bam_file}")
    else:
        logging.info("Mapping completed. Results not saved as SAM file.")
    instream.close()


def cmd_line():
    parser = argparse.ArgumentParser(description="Check PhiX sequences in FASTQ files.")
    parser.add_argument(
        "--R1",
        required=True,
        type=Path,
        help="Path to the R1 FASTQ file.",
    )
    parser.add_argument(
        "--R2",
        required=False,
        type=Path,
        help="Path to the R2 FASTQ file, if applicable.",
    )
    parser.add_argument(
        "--tmpdir",
        required=False,
        type=Path,
        default=gettempdir(),
        help="Temporary directory to store intermediate files like indexes or bam files.",
    )
    parser.add_argument(
        "--outdir",
        required=False,
        type=Path,
        default=Path("."),
        help="Output directory for results.",
    )
    parser.add_argument(
        "--only-estimates",
        required=False,
        action="store_true",
        help="Only estimate the PhiX content without running the full analysis.",
    )
    parser.add_argument(
        "--filter",
        required=False,
        action="store_true",
        help="Filter the reads based on PhiX content. If set, reads that do not map to PhiX will written to the output directory in a file named similarly.",
    )
    parser.add_argument(
        "--keep-alignments",
        required=False,
        action="store_true",
        help="Keep the alignments in SAM format after mapping.",
    )
    parser.add_argument(
        "--loglevel",
        required=False,
        type=str.upper,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level.",
    )
    parser.add_argument(
        "--threads",
        required=False,
        default=1,
        type=int,
        help="Number of threads to use for the mapping process. ",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {distribution('pataphix').version}",
        help="Show the version of pataphix.",
    )
    return parser.parse_args()


def setup_logging(loglevel: str = "INFO"):
    """
    Sets up logging for pataphix.
    """
    str_format = (
        "%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s"  # noqa: E501
    )
    datefmt = "%Y-%m-%d %H:%M:%S"
    logging.basicConfig(
        level=loglevel.upper(),
        format=str_format,
        datefmt=datefmt,
    )
    logging.info(
        f"Command: '{' '.join(sys.argv)}'",
    )
    logging.info("Python version: %s", sys.version)
    logging.info("pataphix version: %s", distribution("pataphix").version)


def process_stderr_log(stderr_output: IO[str] | None, outdir: Path):
    """
    Process the stderr output of the mapping process.
    This function can be extended to handle specific error messages or warnings.
    """
    nb_reads = 0
    perc_reads = "0%"
    all_log = stderr_output.read()
    for line in all_log.splitlines():
        if "overall alignment rate" in line:
            perc_reads = line.split()[0].strip()
            logging.info(f"Overall alignment rate: {perc_reads}")
        if "exactly 1 time" in line:
            nb_reads += int(line.strip().split()[0])
        if ">1 times" in line:
            nb_reads += int(line.strip().split()[0])
    logging.info(f"Number of reads aligned to PhiX: {nb_reads}")
    # write the stderr output to a file
    outdir.mkdir(parents=True, exist_ok=True)
    open(outdir / "mapping_stderr.log", "w").write(all_log)
    fout = open(outdir / "phix.tsv", "w")
    fout.write("phix_reads\tperc_reads\n")
    fout.write(f"{nb_reads}\t{perc_reads}\n")
    fout.close()


def main():
    args = cmd_line()
    tmpdir = TemporaryDirectory(dir=args.tmpdir)
    setup_logging(args.loglevel)
    is_bowtie2_installed()
    index = index_fasta(tmpdir=tmpdir)

    mapping_process = map_reads(
        R1=args.R1,
        R2=args.R2,
        index=index,
        threads=args.threads,
        outdir=args.outdir,
        filter=args.filter,
    )

    if not args.only_estimates:
        logging.info("Processing mapping results...")
        deal_with_mapping(
            process_output=mapping_process.stdout,
            outdir=args.outdir,
            keep=args.keep_alignments,
        )
    mapping_process.wait()
    # process the stderr output of the mapping process
    process_stderr_log(mapping_process.stderr, outdir=args.outdir)
    tmpdir.cleanup()


if __name__ == "__main__":
    main()
