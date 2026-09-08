#!/usr/bin/env python
# Sourced and edited from nf-core/viralrecon:
# https://github.com/nf-core/viralrecon/blob/3731dd3a32a67a2648ea22c2bd980c224abdaee2/bin/multiqc_to_custom_csv.py
import argparse
import csv
import errno
import os
import re
import sys

import yaml


ASSEMBLERS = [
    "autocycler",
    "dragonflye",
    "unicycler",
    "megahit",
    "miniasm",
    "raven",
    "flye",
    "canu",
]

SAMPLE_FIELDS = [
    "Assembly type",
    "# Input short reads",
    "# Trimmed short reads (fastp)",
    "# Input long reads",
    "# Median long read length",
    "# Median long read quality",
    "# Best hit (KmerFinder)",
    "# Best hit assembly ID (KmerFinder)",
    "# Best hit query coverage (KmerFinder)",
    "# Best hit depth (KmerFinder)",
    "# Second hit (KmerFinder)",
    "# Second hit assembly ID (KmerFinder)",
    "# Second hit query coverage (KmerFinder)",
    "# Second hit depth (KmerFinder)",
    "Best assembler",
    "Best assembly # contigs",
    "Best assembly largest contig",
    "Best assembly N50",
    "Best assembly total length",
    "Best assembly GC (%)",
    "Best assembly genome fraction (%)",
]

KMERFINDER_FIELDS = [
    "# Best hit (KmerFinder)",
    "# Best hit assembly ID (KmerFinder)",
    "# Best hit query coverage (KmerFinder)",
    "# Best hit depth (KmerFinder)",
    "# Second hit (KmerFinder)",
    "# Second hit assembly ID (KmerFinder)",
    "# Second hit query coverage (KmerFinder)",
    "# Second hit depth (KmerFinder)",
]

ASSEMBLY_FIELDS = [
    "Sample",
    "Assembly type",
    "Assembler",
    "# Contigs",
    "# Largest contig",
    "# N50",
    "# Total length",
    "# GC (%)",
    "# Genome fraction (%)",
]

SAMPLE_HEADER_CONFIG = {
    "Assembly type": {
        "description": "Read technology inferred for this sample or selected by --assembly_type",
    },
    "# Input short reads": {
        "description": "Total number of short reads in raw FASTQ files",
        "format": "{:,.0f}",
    },
    "# Trimmed short reads (fastp)": {
        "description": "Total number of short reads remaining after adapter and quality trimming with fastp",
        "format": "{:,.0f}",
    },
    "# Input long reads": {
        "description": "Total number of long reads in raw FASTQ files",
        "format": "{:,.0f}",
    },
    "# Median long read length": {
        "description": "Median long-read length in bp",
        "format": "{:,.0f}",
    },
    "# Median long read quality": {
        "description": "Median long-read quality in Phred scale",
        "format": "{:,.1f}",
    },
    "# Best hit (KmerFinder)": {
        "description": "Species name of the best hit from KmerFinder",
    },
    "# Best hit assembly ID (KmerFinder)": {
        "description": "Assembly ID of the best hit from KmerFinder",
    },
    "# Best hit query coverage (KmerFinder)": {
        "description": "Query coverage value of the best hit from KmerFinder",
        "format": "{:,.2f}",
    },
    "# Best hit depth (KmerFinder)": {
        "description": "Depth of the best hit from KmerFinder",
        "format": "{:,.2f}",
    },
    "# Second hit (KmerFinder)": {
        "description": "Species name of the second hit from KmerFinder",
    },
    "# Second hit assembly ID (KmerFinder)": {
        "description": "Assembly ID of the second hit from KmerFinder",
    },
    "# Second hit query coverage (KmerFinder)": {
        "description": "Query coverage value of the second hit from KmerFinder",
        "format": "{:,.2f}",
    },
    "# Second hit depth (KmerFinder)": {
        "description": "Depth of the second hit from KmerFinder",
        "format": "{:,.2f}",
    },
    "Best assembler": {
        "description": "Assembler selected by ranking assemblies from the same sample by QUAST N50. Ties keep the first assembly encountered.",
    },
    "Best assembly # contigs": {
        "description": "Number of contigs for the assembly selected by highest QUAST N50",
        "format": "{:,.0f}",
    },
    "Best assembly largest contig": {
        "description": "Largest contig for the assembly selected by highest QUAST N50",
        "format": "{:,.0f}",
    },
    "Best assembly N50": {
        "description": "Highest assembly N50 observed for this sample",
        "format": "{:,.0f}",
    },
    "Best assembly total length": {
        "description": "Total length for the assembly selected by highest QUAST N50",
        "format": "{:,.0f}",
    },
    "Best assembly GC (%)": {
        "description": "GC percentage for the assembly selected by highest QUAST N50",
        "format": "{:,.2f}",
    },
    "Best assembly genome fraction (%)": {
        "description": "Genome fraction percentage for the assembly selected by highest QUAST N50. This requires QUAST to run with a reference genome.",
        "format": "{:,.2f}",
    },
}

ASSEMBLY_HEADER_CONFIG = {
    "Sample": {
        "description": "Sample identifier",
    },
    "Assembly type": {
        "description": "Read technology inferred for this sample or selected by --assembly_type",
    },
    "Assembler": {
        "description": "Assembler used to generate this assembly",
    },
    "# Contigs": {
        "description": "Total number of contigs calculated by QUAST",
        "format": "{:,.0f}",
    },
    "# Largest contig": {
        "description": "Size of the largest contig calculated by QUAST",
        "format": "{:,.0f}",
    },
    "# N50": {
        "description": "N50 metric for de novo assembly as calculated by QUAST",
        "format": "{:,.0f}",
    },
    "# Total length": {
        "description": "Total assembly length calculated by QUAST",
        "format": "{:,.0f}",
    },
    "# GC (%)": {
        "description": "Assembly GC percentage calculated by QUAST",
        "format": "{:,.2f}",
    },
    "# Genome fraction (%)": {
        "description": "Genome fraction percentage calculated by QUAST. This is only available when QUAST runs with a reference genome.",
        "format": "{:,.2f}",
    },
}


def parse_args(args=None):
    description = "Create custom MultiQC summary tables for nf-core/bacass."
    epilog = "Example usage: multiqc_to_custom_csv.py --assembly_type auto"
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-md",
        "--multiqc_data_dir",
        type=str,
        dest="MULTIQC_DATA_DIR",
        default="multiqc_data",
        help="Directory containing YAML files generated by MultiQC. (default: 'multiqc_data').",
    )
    parser.add_argument(
        "-t",
        "--assembly_type",
        type=str,
        dest="ASSEMBLY_TYPE",
        default="short",
        help="Assembly mode for genome de novo assembly (options: short, long, hybrid, auto).",
    )
    parser.add_argument(
        "-op",
        "--out_prefix",
        type=str,
        dest="OUT_PREFIX",
        default="summary",
        help="Output prefix (default: 'summary').",
    )
    parser.add_argument(
        "-qd",
        "--quast_dir",
        type=str,
        dest="QUAST_DIR",
        default="quast",
        help="Directory containing QUAST report folders staged for MultiQC. (default: 'quast').",
    )
    return parser.parse_args(args)


def make_dir(path):
    if path:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def read_yaml(path):
    if not os.path.exists(path):
        print("WARNING: File does not exist: {}".format(path))
        return {}
    with open(path) as file_handle:
        return yaml.safe_load(file_handle) or {}


def find_tag(data, tag):
    if isinstance(data, dict):
        if tag in data:
            yield data[tag]
        for value in data.values():
            yield from find_tag(value, tag)
    elif isinstance(data, list):
        for value in data:
            yield from find_tag(value, tag)


def first_tag(data, tag):
    values = list(find_tag(data, tag))
    return values[0] if values else None


def nested_value(data, path):
    value = data
    for key in path:
        if not isinstance(value, dict) or key not in value:
            return None
        value = value[key]
    return value


def mapped_value(data, selector):
    if isinstance(selector, (list, tuple)):
        return nested_value(data, selector)
    return first_tag(data, selector)


def kmerfinder_was_run(multiqc_data_dir):
    return os.path.exists(os.path.join(multiqc_data_dir, "multiqc_kmerfinder.yaml"))


def clean_value(value):
    if value is None or value == "":
        return "NA"
    if isinstance(value, str):
        return value.strip() or "NA"
    return value


def numeric_value(value):
    value = clean_value(value)
    if value == "NA":
        return None
    try:
        return float(str(value).replace(",", ""))
    except ValueError:
        return None


def sample_assembly_type(sample, default_assembly_type):
    if default_assembly_type != "auto":
        return default_assembly_type
    if sample.startswith("short-"):
        return "short"
    if sample.startswith("long-"):
        return "long"
    if sample.startswith("hybrid-"):
        return "hybrid"
    return "NA"


def split_assembly_id(assembly_id):
    for assembler in ASSEMBLERS:
        pattern = r"^(?P<sample>.+)-(?P<assembler>{})(?:$|[._-].*)".format(
            re.escape(assembler)
        )
        match = re.match(pattern, assembly_id)
        if match:
            return match.group("sample"), match.group("assembler")
    return assembly_id, "unknown"


def assembly_match_key(assembly_id):
    return split_assembly_id(assembly_id)


def set_metric(rows, key, field, value):
    if key not in rows:
        rows[key] = {}
    rows[key][field] = clean_value(value)


def add_sample_metrics(rows, yaml_dict, mapping):
    for sample, metrics in yaml_dict.items():
        if sample not in rows:
            rows[sample] = {}
        for field, selector in mapping:
            set_metric(rows, sample, field, mapped_value(metrics, selector))


def load_sample_rows(multiqc_data_dir, assembly_type):
    sample_rows = {}

    add_sample_metrics(
        sample_rows,
        read_yaml(os.path.join(multiqc_data_dir, "multiqc_fastp.yaml")),
        [
            ("# Input short reads", ["summary", "before_filtering", "total_reads"]),
            (
                "# Trimmed short reads (fastp)",
                ["summary", "after_filtering", "total_reads"],
            ),
        ],
    )
    add_sample_metrics(
        sample_rows,
        read_yaml(os.path.join(multiqc_data_dir, "multiqc_nanostat.yaml")),
        [
            ("# Input long reads", "Number of reads_fastq"),
            ("# Median long read length", "Median read length_fastq"),
            ("# Median long read quality", "Median read quality_fastq"),
        ],
    )
    add_sample_metrics(
        sample_rows,
        read_yaml(os.path.join(multiqc_data_dir, "multiqc_kmerfinder.yaml")),
        [
            ("# Best hit (KmerFinder)", "07-kmerfinder_best_hit_Species"),
            ("# Best hit assembly ID (KmerFinder)", "07-kmerfinder_best_hit_# Assembly"),
            (
                "# Best hit query coverage (KmerFinder)",
                "07-kmerfinder_best_hit_Query_Coverage",
            ),
            ("# Best hit depth (KmerFinder)", "07-kmerfinder_best_hit_Depth"),
            ("# Second hit (KmerFinder)", "07-kmerfinder_second_hit_Species"),
            (
                "# Second hit assembly ID (KmerFinder)",
                "07-kmerfinder_second_hit_# Assembly",
            ),
            (
                "# Second hit query coverage (KmerFinder)",
                "07-kmerfinder_second_hit_Query_Coverage",
            ),
            ("# Second hit depth (KmerFinder)", "07-kmerfinder_second_hit_Depth"),
        ],
    )

    for sample in sample_rows:
        sample_rows[sample]["Assembly type"] = sample_assembly_type(sample, assembly_type)

    return sample_rows


def load_assembly_rows(multiqc_data_dir, assembly_type):
    assembly_rows = {}
    quast_dict = read_yaml(os.path.join(multiqc_data_dir, "multiqc_quast.yaml"))

    for assembly_id, metrics in quast_dict.items():
        sample, assembler = split_assembly_id(assembly_id)
        assembly_rows[assembly_id] = {
            "Sample": sample,
            "Assembly type": sample_assembly_type(sample, assembly_type),
            "Assembler": assembler,
            "# Contigs": clean_value(first_tag(metrics, "# contigs")),
            "# Largest contig": clean_value(first_tag(metrics, "Largest contig")),
            "# N50": clean_value(first_tag(metrics, "N50")),
            "# Total length": clean_value(first_tag(metrics, "Total length")),
            "# GC (%)": clean_value(first_tag(metrics, "GC (%)")),
            "# Genome fraction (%)": clean_value(first_tag(metrics, "Genome fraction (%)")),
        }

    return assembly_rows


def read_quast_report_tsv(report_tsv):
    with open(report_tsv, newline="") as file_handle:
        reader = csv.reader(file_handle, delimiter="\t")
        rows = list(reader)
    if not rows or rows[0][0] != "Assembly":
        return {}

    assemblies = rows[0][1:]
    metrics_by_assembly = {assembly: {} for assembly in assemblies}
    for row in rows[1:]:
        if not row:
            continue
        metric = row[0]
        for assembly, value in zip(assemblies, row[1:]):
            metrics_by_assembly[assembly][metric] = clean_value(value)
    return metrics_by_assembly


def find_quast_reports(quast_dir):
    if not quast_dir or not os.path.isdir(quast_dir):
        return []
    report_paths = []
    for root, _dirs, files in os.walk(quast_dir, followlinks=True):
        if "report.tsv" in files:
            report_paths.append(os.path.join(root, "report.tsv"))
    return sorted(report_paths)


def add_reference_quast_metrics(assembly_rows, quast_dir):
    rows_by_sample_assembler = {
        assembly_match_key(assembly_id): row for assembly_id, row in assembly_rows.items()
    }

    for report_tsv in find_quast_reports(quast_dir):
        for assembly_id, metrics in read_quast_report_tsv(report_tsv).items():
            genome_fraction = metrics.get("Genome fraction (%)")
            if clean_value(genome_fraction) == "NA":
                continue
            sample, assembler = split_assembly_id(assembly_id)
            row = rows_by_sample_assembler.get((sample, assembler))
            if row:
                row["# Genome fraction (%)"] = clean_value(genome_fraction)


def add_best_assembly_metrics(sample_rows, assembly_rows):
    best_by_sample = {}
    for _assembly_id, row in assembly_rows.items():
        sample = row["Sample"]
        n50 = numeric_value(row.get("# N50"))
        if n50 is None:
            continue
        current = best_by_sample.get(sample)
        current_n50 = numeric_value(current.get("# N50")) if current else None
        if current is None or current_n50 is None or n50 > current_n50:
            best_by_sample[sample] = row

    for sample, assembly_row in best_by_sample.items():
        if sample not in sample_rows:
            sample_rows[sample] = {
                "Assembly type": sample_assembly_type(sample, assembly_row["Assembly type"])
            }
        sample_rows[sample]["Best assembler"] = assembly_row["Assembler"]
        sample_rows[sample]["Best assembly # contigs"] = assembly_row["# Contigs"]
        sample_rows[sample]["Best assembly largest contig"] = assembly_row["# Largest contig"]
        sample_rows[sample]["Best assembly N50"] = assembly_row["# N50"]
        sample_rows[sample]["Best assembly total length"] = assembly_row["# Total length"]
        sample_rows[sample]["Best assembly GC (%)"] = assembly_row["# GC (%)"]
        sample_rows[sample]["Best assembly genome fraction (%)"] = assembly_row[
            "# Genome fraction (%)"
        ]


def table_fields(include_kmerfinder_metrics):
    sample_fields = list(SAMPLE_FIELDS)
    assembly_fields = list(ASSEMBLY_FIELDS)
    if not include_kmerfinder_metrics:
        sample_fields = [field for field in sample_fields if field not in KMERFINDER_FIELDS]
        sample_fields.remove("Best assembly genome fraction (%)")
        assembly_fields.remove("# Genome fraction (%)")
    return sample_fields, assembly_fields


def write_csv(rows, fields, key_header, out_file):
    if not rows:
        return
    make_dir(os.path.dirname(out_file))
    with open(out_file, "w", newline="") as file_handle:
        writer = csv.writer(file_handle)
        writer.writerow([key_header] + fields)
        for key in sorted(rows):
            writer.writerow(
                [key]
                + [rows[key].get(field, "NA") if rows[key].get(field, "NA") != "" else "NA" for field in fields]
            )


def write_custom_yaml(rows, fields, header_config, out_file, section_id, section_name, description):
    if not rows:
        return
    custom_content = {
        "id": section_id,
        "section_name": section_name,
        "description": description,
        "plot_type": "table",
        "data": rows,
        "headers": {field: header_config[field] for field in fields},
    }
    with open(out_file, "w") as file_handle:
        yaml.safe_dump(custom_content, file_handle, sort_keys=False)


def write_table(
    rows,
    fields,
    key_header,
    csv_file,
    yaml_file,
    section_id,
    section_name,
    description,
    header_config,
):
    write_csv(rows, fields, key_header, csv_file)
    write_custom_yaml(
        rows,
        fields,
        header_config,
        yaml_file,
        section_id,
        section_name,
        description,
    )


def main(args=None):
    args = parse_args(args)
    if args.ASSEMBLY_TYPE not in ["short", "long", "hybrid", "auto"]:
        raise ValueError(
            "Unsupported assembly_type '{}'. Expected one of: short, long, hybrid, auto.".format(
                args.ASSEMBLY_TYPE
            )
        )

    sample_rows = load_sample_rows(args.MULTIQC_DATA_DIR, args.ASSEMBLY_TYPE)
    assembly_rows = load_assembly_rows(args.MULTIQC_DATA_DIR, args.ASSEMBLY_TYPE)
    include_reference_metrics = kmerfinder_was_run(args.MULTIQC_DATA_DIR)
    if include_reference_metrics:
        add_reference_quast_metrics(assembly_rows, args.QUAST_DIR)
    add_best_assembly_metrics(sample_rows, assembly_rows)
    sample_fields, assembly_fields = table_fields(include_reference_metrics)

    write_table(
        rows=sample_rows,
        fields=sample_fields,
        key_header="Sample",
        csv_file=args.OUT_PREFIX + "_sample_assembly_metrics.csv",
        yaml_file=args.OUT_PREFIX + "_sample_assembly_metrics_mqc.yaml",
        section_id="sample_summary_assembly_metrics",
        section_name="Sample summary",
        description="One row per sample with read QC, optional KmerFinder taxonomy and the best assembly by N50.",
        header_config=SAMPLE_HEADER_CONFIG,
    )
    write_table(
        rows=assembly_rows,
        fields=assembly_fields,
        key_header="Assembly",
        csv_file=args.OUT_PREFIX + "_comparison_assembly_metrics.csv",
        yaml_file=args.OUT_PREFIX + "_comparison_assembly_metrics_mqc.yaml",
        section_id="assembly_comparison_metrics",
        section_name="Assembly comparison",
        description="One row per assembly so multiple assemblers can be compared without duplicating sample-level QC metrics.",
        header_config=ASSEMBLY_HEADER_CONFIG,
    )


if __name__ == "__main__":
    sys.exit(main())
