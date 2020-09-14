
from unittest import mock
import os.path
import tempfile

from .context import fastq_demux
from fastq_demux.demux import FastqDemultiplexer, Demultiplexer
from fastq_demux.parser import FastqFileParser
from fastq_demux.writer import FastqFileWriter


class TestDemultiplexer:

    def test_get_barcode_file_writers(self, fastq_parser, samplesheet_parser, samplesheet_entries):
        with mock.patch("fastq_demux.demux.FastqFileWriter") as writer_mock:
            demuxer = Demultiplexer(
                fastq_file_parser=fastq_parser,
                samplesheet_parser=samplesheet_parser,
                prefix="this-is-a-prefix",
                outdir="this-is-an-outdir",
                unknown_barcode="this-is-the-unknown-barcode-id")
            file_writers = demuxer.get_barcode_file_writers()
            expected_barcodes = ["+".join(entry.split("\t")[1:]) for entry in samplesheet_entries]
            expected_barcodes.append("this-is-the-unknown-barcode-id")
            assert list(file_writers.keys()) == expected_barcodes

    def test_fastq_file_name_from_sample_id(self, fastq_parser, samplesheet_parser):
        prefix = "this-is-a-prefix-"
        outdir = tempfile.gettempdir()
        barcode = "this-is-a+barcode"
        sample_id = "this-is-a-sample-id"

        fastq_parser.is_single_end = True
        demux = Demultiplexer(
            fastq_file_parser=fastq_parser,
            samplesheet_parser=samplesheet_parser,
            prefix=prefix,
            outdir=outdir,
            unknown_barcode="Unknown")

        expected_filename = os.path.join(
            outdir,
            f"{prefix}{sample_id}_{barcode.replace('+', '-')}_R1.fastq.gz"
        )

        assert demux.fastq_file_name_from_sample_id(sample_id, barcode) == [expected_filename]

        fastq_parser.is_single_end = False
        demux = Demultiplexer(
            fastq_file_parser=fastq_parser,
            samplesheet_parser=samplesheet_parser,
            prefix=prefix,
            outdir=outdir,
            unknown_barcode="Unknown")

        assert demux.fastq_file_name_from_sample_id(sample_id, barcode) == [
            expected_filename,
            expected_filename.replace("_R1.", "_R2.")]


class TestFastqDemultiplexer:

    def test_demultiplex(self, fastq_parser, fastq_writer, fastq_records):
        demuxer = FastqDemultiplexer(fastq_parser, fastq_writer)
        demuxer.demultiplex()
        for i, record in enumerate(fastq_records):
            expected_barcode = record[0][0].split(":")[-1] if i > 0 else "Unknown"
            for j, read_record in enumerate(record):
                fastq_writer[expected_barcode][j].write_record.assert_called_once_with(read_record)

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
