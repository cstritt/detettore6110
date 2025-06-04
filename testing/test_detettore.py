import unittest
import subprocess
import os
import tempfile
import shutil

class TestDetettore6110(unittest.TestCase):
    def setUp(self):
        # Shared inputs
        self.reads = "some_reads.fastq.gz"
        self.targets = "../resources/is_targets/IS6110.fasta"
        self.reference = "../resources/reference/MTBC0_v1.1.fasta"
        self.annotation = "../resources/reference/MTBC0v1.1_PGAP_annot.gff"
        self.prefix = "some_reads"

    def test_run_no_reference(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            result = subprocess.run([
                "../detettore6110.py", self.reads,
                "-t", self.targets,
                "-o", temp_dir,
                "-p", self.prefix
            ], capture_output=True, text=True)

            self.assertEqual(result.returncode, 0, msg=f"Script failed: {result.stderr}")

            anchorfile = os.path.join(temp_dir, f"{self.prefix}.anchors.tsv")
            self.assertTrue(os.path.isfile(anchorfile), "Anchor output missing in no-ref run.")

    def test_run_with_reference(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            result = subprocess.run([
                "../detettore6110.py", self.reads,
                "-t", self.targets,
                "-r", self.reference,
                "-a", self.annotation,
                "-o", temp_dir,
                "-p", self.prefix
            ], capture_output=True, text=True)

            self.assertEqual(result.returncode, 0, msg=f"Script failed: {result.stderr}")

            refins_file = os.path.join(temp_dir, f"{self.prefix}.reference_insertions.tsv")
            anchorfile = os.path.join(temp_dir, f"{self.prefix}.anchors.tsv")  # Adjust if different

            self.assertTrue(os.path.isfile(refins_file), "Reference insertions missing in ref run.")
            self.assertTrue(os.path.isfile(anchorfile), "Anchor output missing in ref run.")

if __name__ == '__main__':
    unittest.main()
