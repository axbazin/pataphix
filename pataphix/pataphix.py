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
from typing import List, IO, Tuple
import os
from collections import defaultdict, Counter

# third-party libraries
import pysam

REFERENCE_FASTA = Path(os.path.dirname(os.path.realpath(__file__)) + "/data/phix.fasta")


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
    bowtie2_cmd: List[str] = [
        "bowtie2-build",
        str(REFERENCE_FASTA),
        tmpdir.name + "/PhiX",
    ]
    logging.debug(f"Launching: {' '.join(bowtie2_cmd)}")
    subprocess.run(bowtie2_cmd, capture_output=True)
    logging.debug(f"Done indexing PhiX fasta file: {REFERENCE_FASTA}")
    logging.debug(f"Index files are in: {tmpdir.name}")
    logging.debug(f"Index files are: {os.listdir(tmpdir.name)}")
    return tmpdir.name + "/PhiX"


def read_fasta(fasta_file: Path) -> str:
    """
    Read the FASTA file and return the sequence as a string.
    This function assumes that the FASTA file is well-formed and contains a single sequence.
    """
    sequence = ""
    with open(fasta_file, "r") as f:
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip()
    return sequence


def get_fasta_length(fasta_file: Path) -> int:
    """
    Get the length of the FASTA file.
    This function assumes that the FASTA file is well-formed and contains a single sequence.
    """
    return len(read_fasta(fasta_file))


def map_reads(
    R1: Path,
    R2: Path | None,
    index: str,
    sam_output: str,
    threads: int,
    outdir: Path,
    filter: bool = False,
) -> subprocess.CompletedProcess[str]:
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
            # This should be improved for cases where '_R1_' is not in the filename
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

    bowtie2_cmd.extend(["-S", sam_output])

    logging.debug(f"Launching: {' '.join(bowtie2_cmd)}")
    mapping_process = subprocess.run(
        bowtie2_cmd,
        check=True,
        capture_output=True,
        text=True,
    )
    return mapping_process


def compute_position_error_rates(
    position_errors: defaultdict[int, int],
    position_coverage: defaultdict[int, int],
    qc_errors: defaultdict[int, int],
    qc_coverage: defaultdict[int, int],
) -> dict[int, dict[str, float]]:
    """
    Calculate position-wise error rates based on the number of errors and coverage.
    This function assumes that position_errors and position_coverage are defaultdicts
    where keys are positions and values are counts of errors or coverage.
    """
    if not position_errors or not position_coverage:
        logging.warning("No position errors or coverage data available.")
        return {}
    # Calculate position-wise error rates
    position_error_rates = {}

    max_position: int = max(position_coverage.keys()) if position_coverage else 0

    for pos in range(max_position + 1):
        if position_coverage[pos] > 0:
            error_rate = position_errors[pos] / position_coverage[pos]
            position_error_rates[pos] = {
                "error_rate": error_rate,
                "errors": position_errors[pos],
                "avg_qc_errors": qc_errors[pos] / position_errors[pos],
                "avg_qc": qc_coverage[pos] / position_coverage[pos],
                "reads": position_coverage[pos],
            }
    return position_error_rates


def process_mapping(
    alignment_file: str,
) -> Tuple[dict[int, Counter], dict[int, Counter]]:
    """
    Process the output of the bowtie2 mapping.
    If keep is True, the SAM file will be read; otherwise, the process output will be read.
    """
    SNV_dict = defaultdict(Counter)
    insertion_dict = defaultdict(Counter)

    instream = pysam.AlignmentFile(alignment_file, "r")
    for read in instream:
        init_insertion = False
        insertion_sequence = ""
        prev_ref_pos = None  # To track the previous reference position for insertions
        # do stuff with the reads
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        # Get aligned pairs: (query_pos, ref_pos, ref_base)
        aligned_pairs = read.get_aligned_pairs()
        for query_pos, ref_pos in aligned_pairs:
            if query_pos is not None and ref_pos is not None:
                SNV_dict[ref_pos][
                    read.query_sequence[query_pos].upper()
                ] += 1  # Count the base at this position

            elif query_pos is None and ref_pos is not None:
                # This is a deletion in the query sequence
                SNV_dict[ref_pos]["-"] += 1
            elif query_pos is not None and ref_pos is None:
                # find the previous ref_pos
                # This is an insertion in the query sequence
                if init_insertion:
                    insertion_sequence += read.query_sequence[query_pos]
                else:
                    init_insertion = True
                    insertion_sequence = read.query_sequence[query_pos]

            if ref_pos is not None:
                if insertion_sequence != "":
                    # reading the insertion is over
                    insertion_dict[prev_ref_pos][insertion_sequence] += 1
                    init_insertion = False
                    insertion_sequence = ""
                prev_ref_pos = ref_pos

    instream.close()
    return SNV_dict, insertion_dict


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
        help="Filter the reads based on PhiX content. If set, reads that do not map to PhiX will be written to the output directory in a file named similarly as the input file(s).",
    )
    parser.add_argument(
        "--keep-alignments",
        required=False,
        action="store_true",
        help="Keep the alignments on the PhiX genome in SAM format after mapping.",
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
        help="Number of threads to use for the mapping process.",
    )
    parser.add_argument(
        "--skip",
        required=False,
        type=int,
        default=10,
        help="Skip the first and last N bases in the reference. Mapping at the end of contigs can be problematic, so this option allows to skip the first and last N bases in the reference.",
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


def process_stderr_log(stderr_output: IO[str], outdir: Path):
    """
    Process the stderr output of the mapping process.
    This function can be extended to handle specific error messages or warnings.
    """
    nb_reads = 0
    perc_reads = "0%"
    for line in stderr_output.splitlines():
        if "overall alignment rate" in line:
            perc_reads = line.split()[0].strip()
            logging.info(f"Overall alignment rate: {perc_reads}")
        if "1 time" in line:
            nb_reads += int(line.strip().split()[0])
        if ">1 times" in line:
            nb_reads += int(line.strip().split()[0])
    logging.info(f"Number of reads (or pairs of reads) aligned to PhiX: {nb_reads}")
    # write the stderr output to a file
    outdir.mkdir(parents=True, exist_ok=True)
    open(outdir / "mapping_stderr.log", "w").write(stderr_output)
    fout = open(outdir / "phix.tsv", "w")
    fout.write("phix_reads\tperc_reads\n")
    fout.write(f"{nb_reads}\t{perc_reads}\n")
    fout.close()


def write_results(results: dict, outdir: Path):
    """
    Write the results of the mapping process to a file in the output directory.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "phix_metrics.tsv", "w") as fout:
        fout.write("metric\tvalue\n")
        for key, value in results.items():
            if not isinstance(value, dict):
                fout.write(f"{key}\t{value}\n")
                logging.debug(f"{key}: {value}")
    logging.info(f"Resulting metrics written to {outdir / 'phix_metrics.tsv'}")
    with open(outdir / "phix_error_estimates.tsv", "w") as fout:
        fout.write(
            "read\tposition\terror_rate\terrors\tavg_qc_errors\tavg_qc\tcoverage\n"
        )
        for pos, metrics in results.get("R1_position_error_rates", {}).items():
            fout.write(
                f"R1\t{pos}\t{metrics['error_rate']:.6f}\t{metrics['errors']}\t{metrics['avg_qc_errors']}\t{metrics['avg_qc']}\t{metrics['reads']}\n"
            )
            logging.debug(
                f"R1\t{pos}\t{metrics['error_rate']:.6f}\t{metrics['errors']}\t{metrics['avg_qc_errors']}\t{metrics['avg_qc']}\t{metrics['reads']}"
            )
        if "R2_position_error_rates" in results:
            for pos, metrics in results["R2_position_error_rates"].items():
                fout.write(
                    f"R2\t{pos}\t{metrics['error_rate']:.6f}\t{metrics['errors']}\t{metrics['avg_qc_errors']}\t{metrics['avg_qc']}\t{metrics['reads']}\n"
                )
                logging.debug(
                    f"R2\t{pos}\t{metrics['error_rate']:.6f}\t{metrics['errors']}\t{metrics['avg_qc_errors']}\t{metrics['avg_qc']}\t{metrics['reads']}"
                )
    logging.info(f"Error estimates written to {outdir / 'phix_error_estimates.tsv'}")


def generate_consensus(
    SNV: defaultdict[int, Counter],
    insertion: defaultdict[int, Counter],
) -> Tuple[list, dict]:
    """
    Generates what would be the 'original' consensus sequence based on the SNV and insertion dictionaries.
    This is a very basic consensus generation that simply takes the most common base at each position.
    It does not take into account the quality of the bases or any other factors.
    It will also give out positions of expected insertions if any.
    """
    consensus = list()
    for char in read_fasta(REFERENCE_FASTA):
        consensus.append(char)  # init consensus with empty strings
    insertion_consensus = {}

    for pos, counts in SNV.items():
        # Get the most common base at this position
        most_common_base, _ = counts.most_common(1)[0]
        if most_common_base != consensus[pos]:
            logging.warning(
                f"High SNV frequency at position {pos}: '{most_common_base}' ({counts[most_common_base]}/{SNV[pos].total()} observations). Using '{most_common_base}' as consensus base instead of '{consensus[pos]}'."
            )
            consensus[pos] = most_common_base

    for pos, ins_counts in insertion.items():
        # Get the most common insertion sequence at this position
        most_common_ins, most_common_count = ins_counts.most_common(1)[0]
        pos_cov = SNV[pos].total()
        if most_common_count > pos_cov * 0.80:
            logging.warning(
                f"High insertion frequency at position {pos}: '{most_common_ins}' ({most_common_count} times, coverage {pos_cov})"
            )
            insertion_consensus[pos] = most_common_ins

    return consensus, insertion_consensus


def compute_error_rates(
    alignment_file: str,
    consensus: list,
    insertion: dict,
    skip: int = 10,
) -> dict:
    """
    Compute the error rates based on the alignment file and the consensus sequence.
    This function calculates the position-wise error rates for both R1 and R2 reads.
    """
    position_errors = {"R1": defaultdict(int), "R2": defaultdict(int)}
    qc_errors = {"R1": defaultdict(int), "R2": defaultdict(int)}
    position_coverage = {"R1": defaultdict(int), "R2": defaultdict(int)}
    qc_coverage = {"R1": defaultdict(int), "R2": defaultdict(int)}

    nb_reads = 0
    nb_single_reads = 0
    nb_paired_reads = 0
    nb_mismatches = 0
    nb_indels = 0

    ref_length = get_fasta_length(REFERENCE_FASTA)

    instream = pysam.AlignmentFile(alignment_file, "r")
    for read in instream:
        nb_reads += 1
        if read.is_paired:
            nb_paired_reads += 1
        else:
            nb_single_reads += 1
        read_type = "R1" if read.is_read1 else "R2"
        init_insertion = False
        insertion_sequence = ""
        prev_ref_pos = None  # To track the previous reference position for insertions
        prev_query_pos = None  # To track the previous query position for insertions
        # do stuff with the reads
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        aligned_pairs = read.get_aligned_pairs()
        for query_pos, ref_pos in aligned_pairs:
            if ref_pos is not None and ref_pos > skip and ref_pos < ref_length - skip:
                # Only consider positions within the reference length and after skipping 'skip' bases
                if query_pos is not None:
                    position_coverage[read_type][query_pos] += 1
                    qc_coverage[read_type][query_pos] += read.query_qualities[query_pos]
                    if read.query_sequence[query_pos] != consensus[ref_pos]:
                        position_errors[read_type][query_pos] += 1
                        qc_errors[read_type][query_pos] += read.query_qualities[
                            query_pos
                        ]
                        nb_mismatches += 1
                else:
                    # there is a deletion in the query sequence
                    # counting the error in the previous position
                    position_coverage[read_type][prev_query_pos] += 1
                    position_errors[read_type][prev_query_pos] += 1
                    nb_indels += 1
                    qc_coverage[read_type][prev_query_pos] += read.query_qualities[
                        prev_query_pos
                    ]
                    qc_errors[read_type][prev_query_pos] += read.query_qualities[
                        prev_query_pos
                    ]
            elif query_pos is not None and ref_pos is None:
                # find the previous ref_pos
                # This is an insertion in the query sequence
                if init_insertion:
                    insertion_sequence += read.query_sequence[query_pos]
                else:
                    init_insertion = True
                    insertion_sequence = read.query_sequence[query_pos]

            if ref_pos is not None:
                if (
                    insertion_sequence != ""
                    and ref_pos > skip
                    and ref_pos < ref_length - skip
                ):
                    # if there is an insertion sequence and we are still in the reference range
                    # then we need to check if the insertion sequence is already in the insertion dict
                    if (not prev_ref_pos in insertion) or (
                        prev_ref_pos in insertion
                        and not insertion_sequence in insertion[prev_ref_pos]
                    ):
                        # then the whole insertion sequence is an error
                        for i in range(prev_query_pos, query_pos):
                            position_errors[read_type][i] += 1
                            nb_indels += 1
                            qc_errors[read_type][i] += read.query_qualities[i]
                    for i in range(prev_query_pos, query_pos):
                        position_coverage[read_type][i] += 1
                        qc_coverage[read_type][i] += read.query_qualities[i]

                    init_insertion = False
                    insertion_sequence = ""
                prev_ref_pos = ref_pos
                if query_pos is not None:
                    prev_query_pos = query_pos

    instream.close()

    # Calculate position-wise error rates
    R1_position_error_rates = compute_position_error_rates(
        position_errors["R1"],
        position_coverage["R1"],
        qc_errors["R1"],
        qc_coverage["R1"],
    )

    result_dict = {
        "R1_position_error_rates": R1_position_error_rates,
        "nb_reads": nb_reads,
        "nb_single_reads": nb_single_reads,
        "nb_paired_reads": nb_paired_reads / 2,  # since paired reads are counted twice
        "nb_mismatches": nb_mismatches,
        "nb_indels": nb_indels,
        "nb_errors": nb_mismatches + nb_indels,
        "nb_errors_R1": sum(
            R1_position_error_rates[pos]["errors"] for pos in R1_position_error_rates
        ),
    }

    if "R2" in position_errors:
        R2_position_error_rates = compute_position_error_rates(
            position_errors["R2"],
            position_coverage["R2"],
            qc_errors["R2"],
            qc_coverage["R2"],
        )
        result_dict["R2_position_error_rates"] = R2_position_error_rates
        result_dict["nb_errors_R2"] = sum(
            R2_position_error_rates[pos]["errors"] for pos in R2_position_error_rates
        )

    return result_dict


def main():
    args = cmd_line()
    tmpdir = TemporaryDirectory(dir=args.tmpdir)
    setup_logging(args.loglevel)
    is_bowtie2_installed()
    index = index_fasta(tmpdir=tmpdir)

    # prepare the output directory
    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.only_estimates:
        sam_output = "/dev/null"
    else:
        # write it in the tmpdir
        sam_output = str(tmpdir.name + "/mapped_phix.sam")
    # Launch the mapping process
    logging.info("Starting the mapping process...")
    mapping_results = map_reads(
        R1=args.R1,
        R2=args.R2,
        sam_output=sam_output,
        index=index,
        threads=args.threads,
        outdir=args.outdir,
        filter=args.filter,
    )
    logging.info("Mapping process completed.")
    if mapping_results.stdout is None:
        raise Exception("Mapping process did not produce any output.")
    if not args.only_estimates:
        if args.keep_alignments:
            logging.info("Sorting and indexing the mapping results...")
            alignment_file = str(args.outdir / "mapped_phix_sorted.bam")
            pysam.sort(
                "--write-index",
                "--threads",
                str(args.threads),
                "-o",
                str(alignment_file),
                str(sam_output),
            )
            logging.info("Mapping results sorted and indexed.")
        else:
            alignment_file = sam_output

        logging.info("Processing the mapping results...")
        # Now deal with the mapping results
        SNV_dict, INS_dict = process_mapping(
            alignment_file=alignment_file,
        )
        logging.info("Mapping results processed. Generating consensus sequence...")
        consensus, insertion_consensus = generate_consensus(
            SNV=SNV_dict,
            insertion=INS_dict,
        )
        logging.info("Consensus sequence generated.")
        logging.info("Calculating read position error rates...")
        results = compute_error_rates(
            alignment_file=alignment_file,
            consensus=consensus,
            insertion=insertion_consensus,
        )
        write_results(results, args.outdir)
        logging.info("Mapping results processed successfully.")

    # process the stderr output of the mapping process
    process_stderr_log(mapping_results.stderr, outdir=args.outdir)
    tmpdir.cleanup()


if __name__ == "__main__":
    main()
