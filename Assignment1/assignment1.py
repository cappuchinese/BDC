"""
Deliverable for first assignment
"""

__author__ = "Lisa Hu"
__date__ = 2023.05


import csv
import multiprocessing as mp
import argparse as ap


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
    if outputcsv is None:
        for i, score in enumerate(scores):
            print(f"{i}: {score}")
    writer = csv.writer(outputcsv)
    for i, score in enumerate(scores):
        writer.writerow([i, score])


def main(argv):
    """

    :param argv:
    :return:
    """
    for fastqfile in argv.fastq_files:
        quality = read_fastq(fastqfile)
        with mp.Pool(processes=argv.n) as pool:
            scores = pool.map(calc_score, quality)
            write_csv(scores, argv.csvfile)


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

    main(args)
