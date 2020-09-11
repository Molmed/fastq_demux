
from unittest import mock
import os.path
import pytest
import tempfile

from .context import fastq_demux
from fastq_demux.demux import FastqDemultiplexer, Demultiplexer
from fastq_demux.parser import FastqFileParser


class TestDemultiplexer:

    def test_fastq_file_name_from_sample_id(self):
        prefix = "this-is-a-prefix-"
        outdir = tempfile.gettempdir()
        barcode = "this-is-a+barcode"
        sample_id = "this-is-a-sample-id"

        demux = Demultiplexer(
            fastq_file_r1="this-is-fastq-r1",
            samplesheet="this-is-the-samplesheet",
            prefix=prefix,
            outdir=outdir,
            unknown_barcode="Unknown")

        expected_filename = os.path.join(
            outdir,
            f"{prefix}{sample_id}_{barcode.replace('+', '-')}_R1.fastq.gz"
        )

        assert demux.fastq_file_name_from_sample_id(sample_id, barcode) == [expected_filename]

        demux = Demultiplexer(
            fastq_file_r1="this-is-fastq-r1",
            fastq_file_r2="this-is-fastq-r2",
            samplesheet="this-is-the-samplesheet",
            prefix=prefix,
            outdir=outdir,
            unknown_barcode="Unknown")

        assert demux.fastq_file_name_from_sample_id(sample_id, barcode) == [
            expected_filename,
            expected_filename.replace("_R1.", "_R2.")]


class TestFastqDemultiplexer:

    def test_demultiplex_record(self, fastq_parser, fastq_writer, fastq_records):

        demuxer = FastqDemultiplexer(fastq_parser, fastq_writer)
        for i, record in enumerate(fastq_records):
            expected_barcode = f"record{i+1}"

            def _write_record(barcode, rec):
                assert barcode == expected_barcode
                assert rec == record

            fastq_writer.write_record = _write_record
            demuxer.demultiplex_record(record)

    def test_demultiplex(self, fastq_parser, fastq_writer, fastq_records):

        def _write_record(barcode, records):
            assert barcode in [record[1][0].split(":")[-1] for record in fastq_records]
            assert records in fastq_records

        fastq_writer.write_record = _write_record
        demuxer = FastqDemultiplexer(fastq_parser, fastq_writer)
        demuxer.demultiplex()

    def test_single_barcode_from_record(self):
        barcode = "ACGTGT"
        record = [
            ":".join(["@", "this", "-is-", "a", " _header_", "w1th[", "s]ome", "0dd", "characters", "and", barcode]),
            "this-is-the-nucleotide-sequence",
            "+",
            "this-is-the-quality-sequence"
        ]
        assert FastqDemultiplexer.barcode_from_record(record) == barcode

    def test_dual_barcode_from_record(self):
        barcode = "ACGTGTAG+TGACATGA"
        record = [
            ":".join(["@", "this", "-is-", "a", " _header_", "w1th[", "s]ome", "0dd", "characters", "and", barcode]),
            "this-is-the-nucleotide-sequence",
            "+",
            "this-is-the-quality-sequence"
        ]
        assert FastqDemultiplexer.barcode_from_record(record) == barcode
