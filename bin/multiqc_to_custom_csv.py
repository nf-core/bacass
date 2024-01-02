#!/usr/bin/env python
# Sourced and Edited from nf-core/viralrecon:
#   https://github.com/nf-core/viralrecon/blob/master/bin/multiqc_to_custom_csv.py#L59
import os
import sys
import errno
import argparse
import yaml


def parse_args(args=None):
    Description = (
        "Create custom spreadsheet for pertinent MultiQC metrics generated by the nf-core/viralrecon pipeline."
    )
    Epilog = "Example usage: python multiqc_to_custom_tsv.py"
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-md",
        "--multiqc_data_dir",
        type=str,
        dest="MULTIQC_DATA_DIR",
        default="multiqc_data",
        help="Full path to directory containing YAML files for each module, as generated by MultiQC. (default: 'multiqc_data').",
    )
    parser.add_argument(
        "-op",
        "--out_prefix",
        type=str,
        dest="OUT_PREFIX",
        default="summary",
        help="Full path to output prefix (default: 'summary').",
    )
    return parser.parse_args(args)


def make_dir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


# Find key in dictionary created from YAML file recursively
# From https://stackoverflow.com/a/37626981
def find_tag(d, tag):
    if tag in d:
        yield d[tag]
    for k, v in d.items():
        if isinstance(v, dict):
            for i in find_tag(v, tag):
                yield i


def yaml_fields_to_dict(yaml_file, append_dict={}, field_mapping_list=[], valid_sample_list=[]):
    integer_fields = [
        "mapped_passed",
        "number_of_SNPs",
        "number_of_indels",
        "MISSENSE",
        "# contigs (>= 0 bp)",
        "# contigs (>= 5000 bp)",
        "Largest contig",
    ]
    if os.path.exists(yaml_file):
        with open(yaml_file) as f:
            yaml_dict = yaml.safe_load(f)
            for k in yaml_dict.keys():
                key = k
                include_sample = True
                if len(valid_sample_list) != 0 and key not in valid_sample_list:
                    include_sample = False
                if include_sample:
                    if key not in append_dict:
                        append_dict[key] = {}
                    if field_mapping_list != []:
                        for i, j in field_mapping_list:
                            val = list(find_tag(yaml_dict[k], j[0]))
                            ## Fix for Cutadapt reporting reads/pairs as separate values
                            if j[0] == "r_written" and len(val) == 0:
                                val = [list(find_tag(yaml_dict[k], "pairs_written"))[0] * 2]
                            if len(val) != 0:
                                val = val[0]
                                if len(j) == 2:
                                    val = list(find_tag(val, j[1]))[0]
                                if j[0] in integer_fields:
                                    val = int(val)
                                if i not in append_dict[key]:
                                    append_dict[key][i] = val
                                else:
                                    print(
                                        "WARNING: {} key already exists in dictionary so will be overwritten. YAML file {}.".format(
                                            i, yaml_file
                                        )
                                    )
                    else:
                        append_dict[key] = yaml_dict[k]
    else:
        print("WARNING: File does not exist: {}".format(yaml_file))
        if len(valid_sample_list) != 0:
            for key in valid_sample_list:
                if key not in append_dict:
                    append_dict[key] = {}
                if field_mapping_list != []:
                    for i, j in field_mapping_list:
                        if i not in append_dict[key]:
                            append_dict[key][i] = "NA"
                        else:
                            print(
                                "WARNING: {} key already exists in dictionary so will be overwritten. YAML file {}.".format(
                                    i, yaml_file
                                )
                            )
                else:
                    append_dict[key] = "NA"
    return append_dict


def metrics_dict_to_file(file_field_list, multiqc_data_dir, out_file, valid_sample_list=[]):
    metrics_dict = {}
    field_list = []
    for yaml_file, mapping_list in file_field_list:
        yaml_file = os.path.join(multiqc_data_dir, yaml_file)
        metrics_dict = yaml_fields_to_dict(
            yaml_file=yaml_file,
            append_dict=metrics_dict,
            field_mapping_list=mapping_list,
            valid_sample_list=valid_sample_list,
        )
        field_list += [x[0] for x in mapping_list]

    if metrics_dict != {}:
        make_dir(os.path.dirname(out_file))
        fout = open(out_file, "w")
        header = ["Sample"] + field_list
        fout.write("{}\n".format(",".join(header)))
        for k in sorted(metrics_dict.keys()):
            row_list = [k]
            for field in field_list:
                if field in metrics_dict[k]:
                    if metrics_dict[k][field]:
                        row_list.append(str(metrics_dict[k][field]).replace(",", ";"))
                    else:
                        row_list.append("NA")
                else:
                    row_list.append("NA")
            fout.write("{}\n".format(",".join(row_list)))
        fout.close()
    return metrics_dict


def main(args=None):
    args = parse_args(args)

    ## File names for MultiQC YAML along with fields to fetch from each file
    illumina_assembly_files = [
        (
            "multiqc_fastp.yaml",
            [
                ("# Input reads", ["before_filtering", "total_reads"]),
                ("# Trimmed reads (fastp)", ["after_filtering", "total_reads"]),
            ]
        ),
        (
            "multiqc_quast_quast_unicycler.yaml",
            [
                ("# Contigs (Unicycler)", ["# contigs (>= 0 bp)"]),
                ("# Largest contig (Unicycler)", ["Largest contig"]),
                ("# N50 (Unicycler)", ["N50"]),
                ("# % Genome fraction", ["Genome fraction (%)"]),
            ],
        ),
        (
            "multiqc_kmerfinder.yaml",
            [
                ("# Best hit (Kmerfinder)", ["07-kmerfinder_best_hit_Species"]),
                ("# Best hit assembly ID (Kmerfinder)", ["07-kmerfinder_best_hit_# Assembly"]),
                ("# Best hit query coverage (Kmerfinder)", ["07-kmerfinder_best_hit_Query_Coverage"]),
                ("# Best hit depth (Kmerfinder)", ["07-kmerfinder_best_hit_Depth"]),
                ("# Second hit (Kmerfinder)", ["07-kmerfinder_second_hit_Species"]),
                ("# Second hit assembly ID (Kmerfinder)", ["07-kmerfinder_second_hit_# Assembly"]),
                ("# Second hit query coverage (Kmerfinder)", ["07-kmerfinder_second_hit_Query_Coverage:"]),
                ("# Second hit depth (Kmerfinder)", ["07-kmerfinder_second_hit_Depth"]),
            ]
        ),
    ]

    ## Write de novo assembly metrics to file
    metrics_dict_to_file(
        file_field_list=illumina_assembly_files,
        multiqc_data_dir=args.MULTIQC_DATA_DIR,
        out_file=args.OUT_PREFIX + "_assembly_metrics_mqc.csv",
        valid_sample_list=[],
    )

if __name__ == "__main__":
    sys.exit(main())
