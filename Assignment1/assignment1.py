"""
Deliverable for first assignment
"""

__author__ = "Lisa Hu"
__date__ = 2023.05


import csv
import multiprocessing as mp
import argparse as ap
import sys


def read_fastq(fastqfile):
    """
    Read the fastq file and calculate the average score for each base position
    :param fastqfile: Path to the fastq file
    :return: List of quality lines
    """
    lines = fastqfile.readlines()
    qual_lines = []

    # Check if valid fastq file
    num_lines = len(lines)
    if num_lines % 4 != 0:
        raise ValueError('Invalid FASTQ format')

    # Get every 4th line (the line with the quality)
    for i in range(3, len(lines), 4):
        qual_lines.append(lines[i].strip())

    return qual_lines


def calc_score(qual_lines):
    """
    Function to calculate the average PHRED score per base position
    :param qual_lines: List of quality lines
    :return:
    """
    # Create a list with a score 0 for all the position
    position_scores = [0] * len(qual_lines[0])
    # Iterate through list of lines
    for line in qual_lines:
        # Enum each character of the line
        for i, char in enumerate(line):
            # Store the sum of the PHRED scores per position
            phred = ord(char) - 33
            position_scores[i] += phred
    # Calculate the average PHRED score per position
    avg_score = [score / len(qual_lines) for score in position_scores]
    return avg_score


def write_csv(scores, outputcsv):
    """
    Write the scores to a CSV file
    :param scores: List of average PHRED scores
    :param outputcsv: Path to output file
    """
    # Check if there was an output file given
    if outputcsv is None:
        # Print the results in stdout
        for i in range(0, len(scores)):
            print(f"{i}: {scores[i]}", file=sys.stdout)
    # Write to csv file when given
    else:
        writer = csv.writer(outputcsv)
        for i in range(0, len(scores)):
            writer.writerow([i, scores[i]])


def main(argv):
    """
    Main function
    :param argv: Command line arguments
    :return:
    """
    if len(argv.fastq_files) > 1:
        # Iterate through all the fastq files given
        for fastqfile in argv.fastq_files:
            # Get the quality lines
            quality = read_fastq(fastqfile)
            # Split filename if in other directory
            if "/" in fastqfile.name:
                filename = fastqfile.name.rsplit("/")[1]
            else:
                filename = fastqfile.name
            # Calculate the scores and write to csv
            with mp.Pool(processes=argv.n) as pool:
                scores = pool.map(calc_score, quality)
                # Write results to output with given name
                with open(f"{filename}.output.csv", "w") as csvfile:
                    write_csv(scores, csvfile)
    else:
        fastqfile = argv.fastq_files[0]
        # Get the quality lines
        quality = read_fastq(fastqfile)
        # Calculate the scores and write to csv
        with mp.Pool(processes=argv.n) as pool:
            scores = pool.map(calc_score, quality)
            write_csv(scores, argv.csvfile)
    return 0


if __name__ == '__main__':
    argparser = ap.ArgumentParser(description="Script voor Opdracht 1 van Big Data Computing")
    argparser.add_argument("-n", action="store",
                           dest="n", required=True, type=int,
                           help="Aantal cores om te gebruiken.")
    argparser.add_argument("-o", action="store", dest="csvfile",
                           type=ap.FileType('w', encoding='UTF-8'),
                           required=False,
                           help="CSV file om de output in op te slaan."
                                "Default is output naar terminal STDOUT")
    argparser.add_argument("fastq_files", action="store", type=ap.FileType('r'), nargs='+',
                           help="Minstens 1 Illumina Fastq Format file om te verwerken")
    args = argparser.parse_args()

    EXITCODE = main(args)
    sys.exit(EXITCODE)
