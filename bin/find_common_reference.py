#!/usr/bin/env python
"""
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII
AUTHOR: Guillermo J. Gorines Cordero
MAIL: guillermo.gorines@urjc.es
VERSION: 0.1
CREATED: Early 2022
REVISED: 18-2-2022
DESCRIPTION:
    Given a directory with kmerfinder results, sum them up
    in an outfile named by the user.

INPUT:
    -DIRECTORY: directory containing all kmerfinder results.
    -OUTFILE: Name of the file to write the whole results in.

OUTPUT:
    -OUTFILE: file containing the kmerfinder results.

USAGE:
    python find_common_reference.py -d [DIRECTORY] -o [OUTFILE]
REQUIREMENTS:
    -Python >= 3.6

DISCLAIMER: This script has been designed for the assembly pipeline of BU-ISCIII.
            Feel free to use it at will, however we dont guarantee its success
            outside its purpose.

================================================================
END_OF_HEADER
================================================================
"""
import os
import sys
import errno
import argparse


def parse_args(args=None):
    """
    Parse the args given to argparser
    """
    Description = "Fetch kmerfinder result files and get the most used reference."
    Epilog = """Example usage: python find_common_reference.py -d <input directory> -o <output file>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("-d", help="Input directory.")
    parser.add_argument("-o", help="Output file.")
    return parser.parse_args(args)


def group_references(kmer_result_dir, out_file):
    """
    Unifies the kmerfinder results, and counts their occurrences
    """
    reference_assembly = {}

    # for file in dir
    for k_file in os.listdir(kmer_result_dir):
        # open file
        with open(os.path.join(kmer_result_dir, k_file), "r") as fh:
            file_lines = fh.readlines()

        # remove heading
        try:
            heading = file_lines[0].split("\t")
            first_line = file_lines[1].split("\t")

            # where is the assembly in the header?
            # find reference according to index
            index_assembly = heading.index("# Assembly")
            reference = first_line[index_assembly]

            # add it to the dict if not there
            if reference not in reference_assembly:
                index_description = heading.index("Description")
                reference_assembly[reference] = [0, first_line[index_description]]
            # sum 1 for another occurrence
            reference_assembly[reference][0] += 1
        except IndexError:
            pass

    # sort it (more occurrences first in file)
    order_reference = dict(
        sorted(reference_assembly.items(), key=lambda x: x[1][0], reverse=True)
    )

    # write it
    with open(out_file, "w") as f_out:
        for key, value in order_reference.items():
            f_out.write(key + "\t" + str(value[0]) + "\t" + value[1] + "\n")
    return


def main(args=None):
    args = parse_args(args)
    group_references(args.d, args.o)


if __name__ == "__main__":
    sys.exit(main())
