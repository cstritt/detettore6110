import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class TestDetettore6110(unittest.TestCase):
    def setUp(self):
        self.repo_root = Path(__file__).resolve().parents[1]
        self.reads = self.repo_root / "testing" / "some_reads.fastq.gz"
        self.targets = self.repo_root / "resources" / "is_targets" / "IS6110.fasta"
        self.reference = self.repo_root / "resources" / "reference" / "MTBC0_v1.1.fasta"
        self.annotation = self.repo_root / "resources" / "reference" / "MTBC0v1.1_PGAP_annot.gff"
        self.prefix = "some_reads"
        self.script = self.repo_root / "detettore6110.py"

    def _run_detettore(self, temp_dir, extra_args):
        command = [
            sys.executable,
            str(self.script),
            str(self.reads),
            "-t",
            str(self.targets),
            "-o",
            temp_dir,
            "-p",
            self.prefix,
        ]
        command.extend(extra_args)
        result = subprocess.run(command, cwd=self.repo_root, capture_output=True, text=True)
        return result

    def test_run_no_reference(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            result = self._run_detettore(temp_dir, [])

            self.assertEqual(result.returncode, 0, msg=f"Script failed: {result.stderr}")

            anchorfile = os.path.join(temp_dir, f"{self.prefix}.anchors.tsv")
            self.assertTrue(os.path.isfile(anchorfile), "Anchor output missing in no-ref run.")

    def test_run_with_reference(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            result = self._run_detettore(temp_dir, [
                "-r",
                str(self.reference),
                "-a",
                str(self.annotation),
            ])

            self.assertEqual(result.returncode, 0, msg=f"Script failed: {result.stderr}")

            refins_file = os.path.join(temp_dir, f"{self.prefix}.reference_insertions.tsv")
            anchorfile = os.path.join(temp_dir, f"{self.prefix}.anchors.tsv")

            self.assertTrue(os.path.isfile(refins_file), "Reference insertions missing in ref run.")
            self.assertTrue(os.path.isfile(anchorfile), "Anchor output missing in ref run.")


if __name__ == "__main__":
    unittest.main()
