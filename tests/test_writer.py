
import os.path
import pytest
import tempfile
from unittest import mock

from .context import fastq_demux
from fastq_demux.writer import FastqFileWriter, FastqFileWriterHandler
from fastq_demux.parser import FastqFileParser
from fastq_demux.demux import FastqDemultiplexer


class TestFastqFileWriter:

    def test_write_record(self, fastq_records, fastq_output_file_r1):
        writer = FastqFileWriter(fastq_output_file_r1)
        for record in fastq_records:
            writer.write_record(record[0])
        writer.close_handle()
        parser = FastqFileParser(fastq_output_file_r1)
        assert list(parser.fastq_records()) == [[record[0]] for record in fastq_records]


class TestFastqFileWriterHandler:

    def test_fastq_file_name_from_sample_id(self):
        prefix = "this-is-a-prefix-"
        outdir = tempfile.gettempdir()
        barcode = "this-is-a+barcode"
        sample_id = "this-is-a-sample-id"

        fastq_file_writer_handler = FastqFileWriterHandler(
            prefix=prefix,
            outdir=outdir,
            is_single_end=True,
            no_gzip_compression=True)

        expected_filename = os.path.join(
            outdir,
            f"{prefix}{sample_id}_{barcode.replace('+', '-')}_R1.fastq"
        )

        assert fastq_file_writer_handler.fastq_file_name_from_sample_id(
            sample_id, barcode) == [expected_filename]

        fastq_file_writer_handler.is_single_end = False
        fastq_file_writer_handler.no_gzip_compression = False

        expected_filename = f"{expected_filename}.gz"

        assert fastq_file_writer_handler.fastq_file_name_from_sample_id(
            sample_id, barcode) == [
                expected_filename,
                expected_filename.replace("_R1.", "_R2.")]

    def test_fastq_file_writers_from_sample_id_se(self):
        handler = FastqFileWriterHandler(
            prefix="", outdir="", is_single_end=True, no_gzip_compression=False)
        self._helper_fastq_file_writers_from_sample_id(handler)

    def test_fastq_file_writers_from_sample_id_pe(self):
        handler = FastqFileWriterHandler(
            prefix="", outdir="", is_single_end=False, no_gzip_compression=False)
        self._helper_fastq_file_writers_from_sample_id(handler)

    def _helper_fastq_file_writers_from_sample_id(self, handler):
        with mock.patch("fastq_demux.writer.FastqFileWriter", autospec=True) as writer_mock:
            sample_id = "this-is-the-sample_id"
            barcode = "this-is-the-barcode"
            handler.fastq_file_writers_from_sample_id(sample_id=sample_id, barcode=barcode)
            assert len(writer_mock.call_args_list) == (2 - int(handler.is_single_end))
            for read_no, call_args in enumerate(writer_mock.call_args_list):
                assert sample_id in call_args[0][0]
                assert barcode in call_args[0][0]
                assert f"_R{read_no+1}" in call_args[0][0]

    def test_fastq_file_writers_from_mapping(self, samplesheet_parser):
        handler = FastqFileWriterHandler(
            prefix="", outdir="", is_single_end=False, no_gzip_compression=False)

        def _fastq_file_writers_from_sample_id(*args):
            return "-".join(args)

        with mock.patch.object(
                handler,
                "fastq_file_writers_from_sample_id",
                _fastq_file_writers_from_sample_id):
            barcode_to_sample_mapping = samplesheet_parser.get_barcode_to_sample_mapping()
            assert handler.fastq_file_writers is None
            writers = handler.fastq_file_writers_from_mapping(barcode_to_sample_mapping)
            assert handler.fastq_file_writers == writers
            for barcode, writer in writers.items():
                assert barcode in barcode_to_sample_mapping
                assert writer == "-".join([barcode_to_sample_mapping[barcode], barcode])

    def test_write_fastq_record(self, fastq_writer, fastq_records):
        # make sure that an unknown barcode throws a KeyError
        with pytest.raises(KeyError):
            barcode = FastqDemultiplexer.barcode_from_record(fastq_records[0][0])
            fastq_writer.write_fastq_record(barcode, fastq_records[0])

        # write the unknown barcode record using the "Unknown" barcode
        fastq_writer.write_fastq_record("Unknown", fastq_records[0])
        assert all(
            [writer.write_record.call_count == 1
             for writer in fastq_writer.fastq_file_writers["Unknown"]])

        # the rest should have known barcodes
        for record in fastq_records[1:]:
            barcode = FastqDemultiplexer.barcode_from_record(record[0])
            fastq_writer.write_fastq_record(barcode, record)
            assert all(
                [writer.write_record.call_count == 1
                 for writer in fastq_writer.fastq_file_writers[barcode]])
