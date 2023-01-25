
import json
import os
import tempfile

from click.testing import CliRunner
import pytest

from .context import fastq_demux
from fastq_demux import fastq_demux


class TestFastqDemux:

    @staticmethod
    def _parse_stats_json(statsfile):
        with open(statsfile) as fh:
            stats = json.load(fh)
            return stats

    @staticmethod
    def _demultiplex_helper(
            expected_samples,
            ssheet_file,
            fq_file_r1,
            fq_file_r2=None,
            fq_file_i1=None,
            fq_file_i2=None,
            mismatches=0):
        prefix = "this-is-a-prefix"
        unknown_barcode = "this-is-the-unknown-prefix"
        total_expected_reads = 50
        expected_samples[unknown_barcode] = [
            "barcode",
            (total_expected_reads - sum([s[1] for s in expected_samples.values()])),
            (100.0 - sum([s[2] for s in expected_samples.values()]))]
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
                f"--mismatches={mismatches} "
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

            statsfile = f"{prefix}demux_Stats.json"
            assert statsfile in os.listdir(outdir)

            observed_samples = {}
            stats = TestFastqDemux._parse_stats_json(os.path.join(outdir, statsfile))
            for sample in stats.get("ConversionResults", {}).get("DemuxResults", {}):
                for index_metric in sample["IndexMetrics"]:
                    barcode = index_metric["IndexSequence"]
                    observed_samples[barcode] = [
                        sample["SampleId"],
                        sum(list(index_metric["MismatchCounts"].values())),
                        0]

            observed_samples[unknown_barcode] = [
                "barcode",
                stats.get("ConversionResults", {}).get("Undetermined", {}).get("NumberReads"),
                0]

            total_observed_reads = sum([s[1] for s in observed_samples.values()])
            for sample in observed_samples.keys():
                observed_samples[sample][2] = \
                    100.0 * int(observed_samples[sample][1])/total_observed_reads

            assert total_observed_reads == total_expected_reads
            assert observed_samples == expected_samples

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

    def test_demultiplex_dual_index_mismatch(
            self,
            dual_index_samplesheet_file,
            fastq_file_r1,
            fastq_file_r2,
            fastq_file_i1,
            fastq_file_i2):
        expected_samples = {
            "GGGGGGGG+AGATCTCG": ["Sample1", 25, 50.0],
            "GAAGATTT+TTTACTCT": ["Sample2", 6, 12.0],
            "GAAGATTT+AAAACGCC": ["Sample3", 3, 6.0]}
        self._demultiplex_helper(
            expected_samples,
            dual_index_samplesheet_file,
            fq_file_r1=fastq_file_r1,
            fq_file_r2=fastq_file_r2,
            fq_file_i1=fastq_file_i1,
            fq_file_i2=fastq_file_i2,
            mismatches=1)

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

