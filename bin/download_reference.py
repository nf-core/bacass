#!/usr/bin/env python
"""
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII
AUTHOR: Guillermo J. Gorines Cordero
EDITED BY: Daniel VM
VERSION: 0.1
CREATED: Early 2022
REVISED: 18-2-2022
EDITED: 14-11-2023
DESCRIPTION: 20-05-2024
    Given a file with the kmerfinder results and frequencies (probably
    created by find_common_reference.py), and the NCBI assembly sheet,
    download the top-reference genome, gff and protein files from
    the NCBI ftp.

INPUT:
    -FILE: file containing the ranking of references from kmerfinder created by the script find_common_references
    -REFERENCE: file with the NCBI reference list
    -OUTDIR: name of the output dir

OUTPUT:
    - *_fna.gz: file with the top-reference genome
    - *_gff.gz: file with the top-reference gff
    - *_protein.gz: file with the top-reference proteins

USAGE:
    python download_reference.py
    -file [FILE]
    -reference [REFERENCE]
    -out_dir [OUTDIR]

REQUIREMENTS:
    -Python >= 3.6
    -Python wget

DISCLAIMER:
    This script has been designed for the assembly pipeline of BU-ISCIII.
    Feel free to use it at will, however we dont guarantee its success
    outside its purpose.
================================================================
END_OF_HEADER
================================================================
"""

import sys
import argparse
import os

# import wget
import requests


# TODO: Generate report
def parse_args(args=None):
    Description = "download the reference files \
        (fna, faa, gff)from the reference NCBI file."
    Epilog = """Usage example: \
        python download_reference.py \
        -file <file with the references created by find_common_reference> \
        -reference <file from the NCBI with all bacterial references> \
        -out_dir <output directory>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-file", help="File containing the ranking of references from kmerfinder."
    )
    parser.add_argument(
        "-reference",
        help="File containing the paths to bacterial references. See example in: https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt",
    )
    parser.add_argument("-out_dir", help="Output directory.")

    return parser.parse_args(args)


def download_references(file, reference, out_dir):
    """
    Downloads the top reference from the NCBI database
    """

    reference_ends = ["_genomic.fna.gz", "_protein.faa.gz", "_genomic.gff.gz"]

    # extract the most common reference from file
    with open(file) as infile:
        infile = infile.readlines()
        infile = [
            item.replace("\n", "").split("\t")
            for item in infile
            if not item.startswith("#")
        ]
        top_reference = infile[0][0]

    with open(str(top_reference) + ".winner", "w") as topref:
        topref.write(top_reference)

    # create the outdir (do nothing if already there)
    try:
        os.mkdir(out_dir)
    except FileExistsError:
        pass

    # open the reference and find the reference
    with open(reference) as inref:
        inref = inref.readlines()
        inref = [
            item.replace("\n", "").split("\t")
            for item in inref
            if not item.startswith("#")
        ]

        # Initialize an empty list to store the URLs
        dir_url = []

        # Iterate over each row in the inref
        for row in inref:
            # Construct the ref_query using assembly_accession and asm_name
            assembly_accession = row[0]
            asm_name = row[15]
            ref_query = f"{assembly_accession}_{asm_name}"

            # Check if ref_query matches the search value
            if ref_query == top_reference:
                # make url  # Append the 20th element of the row to the URL list:
                assembly_url = row[19] + "/" + ref_query
                dir_url.append(assembly_url)

        if len(dir_url) == 0:
            print(
                "No assemblies responding to the top reference: ",
                top_reference,
                " were found",
            )
            sys.exit(1)

        dir_url = str(dir_url[0])

    # get url and reference file
    for r_end in reference_ends:
        out_file = out_dir + "/" + top_reference + r_end
        file_url = dir_url + r_end
        print(file_url)

        # wget.download(file_url, out_file)
        response = requests.get(file_url, stream=True)
        with open(out_file, "wb") as out:
            for chunk in response.iter_content(chunk_size=8192):
                out.write(chunk)

    return


def main(args=None):
    args = parse_args(args)
    download_references(args.file, args.reference, args.out_dir)


if __name__ == "__main__":
    sys.exit(main())
