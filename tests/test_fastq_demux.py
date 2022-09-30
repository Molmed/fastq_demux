
import os
import tempfile

from click.testing import CliRunner
import pytest

from .context import fastq_demux
from fastq_demux import fastq_demux


class TestFastqDemux:

    @staticmethod
    def _demultiplex_helper(
            expected_samples,
            ssheet_file,
            fq_file_r1,
            fq_file_r2=None,
            fq_file_i1=None,
            fq_file_i2=None):
        prefix = "this-is-a-prefix"
        unknown_barcode = "this-is-the-unknown-prefix"
        expected_samples[unknown_barcode] = ["barcode"]
        no_reads = int(fq_file_r1 is not None) + int(fq_file_r2 is not None)
        with tempfile.TemporaryDirectory(prefix="TestFastqDemux") as outdir:
            reads = ""
            for flag, fqfile in zip(("--R2", "--I1", "--I2"), (fq_file_r2, fq_file_i1, fq_file_i2)):
                if fqfile:
                    reads = f"{reads} {flag}={fqfile}"

            args = (
                f"--R1={fq_file_r1} "
                f"{reads} "
                f"--samplesheet={ssheet_file} "
                f"--prefix={prefix} "
                f"--outdir={outdir} "
                f"--unknown-barcode={unknown_barcode} "
                f"--no-gzip-compression"
            )
            runner = CliRunner()
            result = runner.invoke(
                fastq_demux.demultiplex,
                args=args)

            assert result.exit_code == 0

            expected_fastq_files = [
                f"{prefix}{sample[0]}_{barcode.replace('+', '-')}_R{read_no+1}.fastq"
                for read_no in range(no_reads)
                for barcode, sample in expected_samples.items()]
            assert sorted(
                expected_fastq_files) == sorted(
                filter(lambda x: x.endswith(".fastq"), os.listdir(outdir)))

            assert f"{prefix}demux_Stats.json" in os.listdir(outdir)

            for line in result.stdout.split("\n"):
                if not line:
                    continue
                barcode, count, pct = line.split()
                if barcode in expected_samples:
                    assert count == str(expected_samples[barcode][1]) and \
                           pct == str(expected_samples[barcode][2])

    def test_demultiplex_single_index(
            self,
            single_index_samplesheet_file,
            fastq_file_r1,
            fastq_file_r2,
            fastq_file_i1):
        expected_samples = {
            "GGGGGGGG": ["Sample1", 26, 52.0],
            "GAAGATTT": ["Sample2", 12, 24.0]}
        self._demultiplex_helper(
            expected_samples,
            single_index_samplesheet_file,
            fq_file_r1=fastq_file_r1,
            fq_file_r2=fastq_file_r2,
            fq_file_i1=fastq_file_i1)

    def test_demultiplex_dual_index(
            self,
            dual_index_samplesheet_file,
            fastq_file_r1,
            fastq_file_r2,
            fastq_file_i1,
            fastq_file_i2):
        expected_samples = {
            "GGGGGGGG+AGATCTCG": ["Sample1", 21, 42.0],
            "GAAGATTT+TTTACTCT": ["Sample2", 6, 12.0],
            "GAAGATTT+AAAACGCC": ["Sample3", 2, 4.0]}
        self._demultiplex_helper(
            expected_samples,
            dual_index_samplesheet_file,
            fq_file_r1=fastq_file_r1,
            fq_file_r2=fastq_file_r2,
            fq_file_i1=fastq_file_i1,
            fq_file_i2=fastq_file_i2)

    def test_demultiplex_single_index_header(
            self,
            single_index_samplesheet_file,
            single_index_fastq_file_r1):
        expected_samples = {
            "GGGGGGGG": ["Sample1", 27, 54.0],
            "GAAGATTT": ["Sample2", 11, 22.0]}
        self._demultiplex_helper(
            expected_samples,
            single_index_samplesheet_file,
            fq_file_r1=single_index_fastq_file_r1)

    def test_demultiplex_dual_index_header(
            self,
            dual_index_samplesheet_file,
            fastq_file_r1,
            fastq_file_r2):
        expected_samples = {
            "GGGGGGGG+AGATCTCG": ["Sample1", 22, 44.0],
            "GAAGATTT+TTTACTCT": ["Sample2", 5, 10.0],
            "GAAGATTT+AAAACGCC": ["Sample3", 3, 6.0]}
        self._demultiplex_helper(
            expected_samples,
            dual_index_samplesheet_file,
            fq_file_r1=fastq_file_r1,
            fq_file_r2=fastq_file_r2)

    def test_demultiplex_illegal_option(
            self,
            dual_index_samplesheet_file,
            fastq_file_r1,
            fastq_file_r2,
            fastq_file_i1,
            fastq_file_i2):
        with pytest.raises(Exception):
            self._demultiplex_helper(
                {},
                dual_index_samplesheet_file,
                fq_file_r1=fastq_file_r1,
                fq_file_r2=fastq_file_r2,
                fq_file_i2=fastq_file_i2)

