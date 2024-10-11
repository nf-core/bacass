#!/usr/bin/env python
"""
Author: Daniel VM
Email: da.valle@ciberisciii.es
Date: 2024/01/20

MIT License

© 2024 Daniel VM

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
"""

import sys
import argparse
import csv
import yaml


def parse_args(args=None):
    Description = "Create a yaml file from csv input file grouping samples as keys and resting fields as their value pair."

    Epilog = "Example usage: python csv_to_yaml.py -i myfile.csv -k 'sample_name' -o converted_file"
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-i", "--input", type=str, dest="CSV_FILE", help="Input file in CSV format."
    )

    parser.add_argument(
        "-k",
        "--key_field",
        type=str,
        dest="KEY_FIELD",
        help="Name of the key/column grupping field in the input csv.",
    )

    parser.add_argument(
        "-op",
        "--output_prefix",
        type=str,
        default="output_file",
        dest="OUT_PREFIX",
        help="Output file name",
    )
    return parser.parse_args(args)


def parse_csv(csv_file):
    with open(csv_file, "r") as c:
        csv_reader = csv.DictReader(c)
        data = [row for row in csv_reader]
    return data


def create_yaml(data, key, output_prefix):
    yaml_data = {
        entry[key]: {k: v for k, v in entry.items() if k != key} for entry in data
    }
    with open(output_prefix + ".yaml", "w") as yaml_file:
        yaml.dump(yaml_data, yaml_file, default_flow_style=False)


def main(args=None):
    args = parse_args(args)
    file_list = parse_csv(args.CSV_FILE)

    create_yaml(data=file_list, key=args.KEY_FIELD, output_prefix=args.OUT_PREFIX)


if __name__ == "__main__":
    sys.exit(main())
