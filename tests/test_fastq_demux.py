
import os
import tempfile

from click.testing import CliRunner

from .context import fastq_demux
from fastq_demux import fastq_demux


class TestFastqDemux:

    def test_demultiplex(self, fastq_file_r1, fastq_file_r2, samplesheet_file):
        prefix = "this-is-a-prefix"
        unknown_barcode = "this-is-the-unknown-prefix"
        expected_samples = {
            "GGGGGGGG+AGATCTCG": ["Sample1", 22, 73.3],
            "GAAGATTT+TTTACTCT": ["Sample2", 5, 16.7],
            "GAAGATTT+AAAACGCC": ["Sample3", 3, 10.0],
            unknown_barcode: ["Sample"]}

        with tempfile.TemporaryDirectory(prefix="TestFastqDemux") as outdir:
            args = (
                f"--R1={fastq_file_r1} "
                f"--R2={fastq_file_r2} "
                f"--samplesheet={samplesheet_file} "
                f"--prefix={prefix} "
                f"--outdir={outdir} "
                f"--unknown-barcode={unknown_barcode}"
            )
            runner = CliRunner()
            result = runner.invoke(
                fastq_demux.demultiplex,
                args=args)

            assert result.exit_code == 0

            expected_fastq_files = [
                f"{prefix}{sample[0]}_{barcode.replace('+', '-')}_R{read_no+1}.fastq.gz"
                for read_no in range(2)
                for barcode, sample in expected_samples.items()]
            assert sorted(expected_fastq_files) == sorted(os.listdir(outdir))

            for line in result.stdout.split("\n"):
                if not line:
                    continue
                barcode, count, pct = line.split()
                if barcode in expected_samples:
                    assert count == str(expected_samples[barcode][1]) and \
                           pct == str(expected_samples[barcode][2])
