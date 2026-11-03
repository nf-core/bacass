import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import yaml


PROJECT_DIR = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_DIR / "bin" / "multiqc_to_custom_csv.py"


class MultiqcToCustomCsvTest(unittest.TestCase):
    def test_reference_metrics_are_retained_without_kmerfinder(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            multiqc_data_dir = temp_path / "multiqc_data"
            multiqc_data_dir.mkdir()
            quast_dir = temp_path / "quast" / "sample"
            quast_dir.mkdir(parents=True)

            with (multiqc_data_dir / "multiqc_quast.yaml").open("w") as file_handle:
                yaml.safe_dump(
                    {
                        "sample-unicycler": {
                            "# contigs": "10",
                            "Largest contig": "1000",
                            "N50": "800",
                            "Total length": "5000",
                            "GC (%)": "50.0",
                        }
                    },
                    file_handle,
                )

            (quast_dir / "report.tsv").write_text(
                "Assembly\tsample-unicycler\n"
                "# contigs (>= 0 bp)\t10\n"
                "Genome fraction (%)\t97.1\n"
            )

            output_prefix = temp_path / "summary"
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--multiqc_data_dir",
                    str(multiqc_data_dir),
                    "--quast_dir",
                    str(temp_path / "quast"),
                    "--out_prefix",
                    str(output_prefix),
                    "--assembly_type",
                    "auto",
                ],
                check=True,
            )

            with (temp_path / "summary_comparison_assembly_metrics.csv").open() as file_handle:
                assembly_rows = list(csv.DictReader(file_handle))
            with (temp_path / "summary_sample_assembly_metrics.csv").open() as file_handle:
                sample_rows = list(csv.DictReader(file_handle))

            self.assertEqual(assembly_rows[0]["# Genome fraction (%)"], "97.1")
            self.assertEqual(sample_rows[0]["Best assembly genome fraction (%)"], "97.1")
            self.assertNotIn("# Best hit (KmerFinder)", sample_rows[0])


if __name__ == "__main__":
    unittest.main()
