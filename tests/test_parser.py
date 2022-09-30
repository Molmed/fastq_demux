
import gzip
import pytest

from .context import fastq_demux
from fastq_demux.parser import FastqFileParser, FastqFileParserHeaderIndex, FastqFileParserI1, \
    FastqFileParserI2, SampleSheetParser
from io import StringIO
from unittest import mock


class TestFastqFileParser:

    def test_factory_method(self):
        assert isinstance(
            FastqFileParser.create_fastq_file_parser(fastq_r1="r1"), FastqFileParserHeaderIndex)
        assert isinstance(
            FastqFileParser.create_fastq_file_parser(fastq_r1="r1", fastq_r2="r2"),
            FastqFileParserHeaderIndex)
        assert isinstance(
            FastqFileParser.create_fastq_file_parser(fastq_r1="r1", fastq_i1="i1"),
            FastqFileParserI1)
        assert isinstance(
            FastqFileParser.create_fastq_file_parser(fastq_r1="r1", fastq_i1="i1", fastq_i2="i2"),
            FastqFileParserI2)

    def test_parser(self, fastq_file_r1, fastq_file_r2):

        # just iterate over the example file and assert that the correct number of records were parsed
        parser = FastqFileParser.create_fastq_file_parser(fastq_file_r1, fastq_file_r2)
        records = list(parser.fastq_records())
        assert len(records) == 50

    def test_parser_single_end(self, fastq_file_r1):

        # iterate over the file and make sure nothing is broken
        parser = FastqFileParser.create_fastq_file_parser(fastq_file_r1)
        assert len(list(parser.fastq_records())) == 50

    def test_parser_mismatching_files(self, fastq_file_r1, fastq_file_short_r2):
        # iterate over the mismatching files and make sure an exception is thrown
        parser = FastqFileParser.create_fastq_file_parser(fastq_file_r1, fastq_file_short_r2)
        with pytest.raises(Exception):
            list(parser.fastq_records())

    def test_ensure_filehandles_empty(self):

        # assert that a single-end parser returns true
        parser = FastqFileParser("this-is-fastq-r1")
        handles = [StringIO(), StringIO()]
        assert parser.ensure_filehandles_empty(handles)

        # assert that a paired-end parser returns true for empty handles
        parser = FastqFileParser("this-is-fastq-r1", "this-is-fastq-r2")
        assert parser.ensure_filehandles_empty(handles)

        # assert that a paired-end parser returns false for non-empty handle
        assert not parser.ensure_filehandles_empty([
            StringIO("this.is.some.contents\n"),
            StringIO()])
        # assert that a paired-end parser returns false for non-empty handle
        assert not parser.ensure_filehandles_empty([
            StringIO(),
            StringIO("this.is.some.contents\n")])

    def test_open_func(self):
        input_output = {
            "this-is-a-gzipped-file.gz": gzip.open,
            "this-is-another-gzipped-file.gzip": gzip.open,
            "this-is-a-plain-file": open
        }
        for input_file, expected_func in input_output.items():
            assert FastqFileParser.open_func(input_file) == expected_func

    def test__file_parser_handle(self, fastq_file_r1):
        parser = FastqFileParser.create_fastq_file_parser(fastq_r1=fastq_file_r1)
        line = next(parser._file_parser_handle(parser.fastq_r1))
        assert line.startswith("@") and len(line.split(":")) == 10


class TestFastqFileParserHeaderIndex:

    def test_single_barcode_from_record(self, fastq_record_single_index_header):
        assert FastqFileParserHeaderIndex.barcode_from_record(
            [fastq_record_single_index_header[0]]) == fastq_record_single_index_header[1]

    def test_dual_barcode_from_record(self, fastq_record_dual_index_header):
        assert FastqFileParserHeaderIndex.barcode_from_record(
            [fastq_record_dual_index_header[0]]) == fastq_record_dual_index_header[1]

    def test_records_from_handles(
            self,
            paired_end_single_index_fastq_records,
            paired_end_single_index_fastq_handles):
        parser = FastqFileParser.create_fastq_file_parser(fastq_r1="r1", fastq_r2="r2")
        observed_fastq_records = []
        while True:
            try:
                observed_fastq_records.append(
                    parser.records_from_handles(
                        paired_end_single_index_fastq_handles[0],
                        paired_end_single_index_fastq_handles[1]))
            except StopIteration:
                break
        assert observed_fastq_records == list(paired_end_single_index_fastq_records)


class TestFastqFileParserI1:

    def test_barcode_from_record(self, fastq_record_single_index_header):
        assert FastqFileParserI1.barcode_from_record(
            [fastq_record_single_index_header[0]]) == fastq_record_single_index_header[0][1]

    def test_records_from_handles(
            self,
            paired_end_single_index_fastq_records,
            paired_end_single_index_fastq_handles):
        parser = FastqFileParser.create_fastq_file_parser(
            fastq_r1="r1", fastq_r2="r2", fastq_i1="i1")
        observed_records = []
        while True:
            try:
                observed_records.append(
                    parser.records_from_handles(
                        paired_end_single_index_fastq_handles[0],
                        paired_end_single_index_fastq_handles[1]))
            except StopIteration:
                break
        assert observed_records == list(paired_end_single_index_fastq_records)


class TestFastqFileParserI2:

    def test_barcode_from_record(
            self,
            fastq_record_single_index_header,
            fastq_record_dual_index_header):
        assert FastqFileParserI2.barcode_from_record(
            [fastq_record_single_index_header[0], fastq_record_dual_index_header[0]]) == \
               "+".join([
                   fastq_record_single_index_header[0][1],
                   fastq_record_dual_index_header[0][1]])

    def test_records_from_handles(
            self,
            paired_end_dual_index_fastq_records,
            paired_end_dual_index_fastq_handles):
        parser = FastqFileParser.create_fastq_file_parser(
            fastq_r1="r1", fastq_r2="r2", fastq_i1="i1", fastq_i2="i2")
        observed_records = []
        while True:
            try:
                observed_records.append(
                    parser.records_from_handles(
                        paired_end_dual_index_fastq_handles[0],
                        paired_end_dual_index_fastq_handles[1]))
            except StopIteration:
                break
        assert observed_records == list(paired_end_dual_index_fastq_records)


class TestSampleSheetParser:

    def test_get_barcode_to_sample_mapping(self, samplesheet_entries):
        parser = SampleSheetParser("this-is-a-samplesheet")
        with mock.patch.object(parser, "get_file_handle", side_effect=[samplesheet_entries]):
            mapping = list()
            for barcode, sample_id in parser.get_barcode_to_sample_mapping().items():
                mapping.append("\t".join([sample_id] + barcode.split("+")))
            assert sorted(mapping) == sorted(samplesheet_entries)

    def test_samplesheet_parser(self, dual_index_samplesheet_file):
        parser = SampleSheetParser(dual_index_samplesheet_file)
        mapping = parser.get_barcode_to_sample_mapping()
        expected_mapping = {
            "GGGGGGGG+AGATCTCG": "Sample1",
            "GAAGATTT+TTTACTCT": "Sample2",
            "GAAGATTT+AAAACGCC": "Sample3"}
        assert mapping == expected_mapping
