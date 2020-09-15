
import gzip
import pytest

from .context import fastq_demux
from fastq_demux.parser import FastqFileParser, SampleSheetParser
from io import StringIO
from unittest import mock


class TestFastqFileParser:

    def test_parser(self, fastq_file_r1, fastq_file_r2):

        # just iterate over the example file and assert that the correct number of records were parsed
        parser = FastqFileParser(fastq_file_r1, fastq_file_r2)
        assert len(list(parser.fastq_records())) == 50

    def test_parser_single_end(self, fastq_file_r1):

        # iterate over the file and make sure nothing is broken
        parser = FastqFileParser(fastq_file_r1)
        assert len(list(parser.fastq_records())) == 50

    def test_parser_mismatching_files(self, fastq_file_r1, fastq_file_short_r2):
        # iterate over the mismatching files and make sure an exception is thrown
        parser = FastqFileParser(fastq_file_r1, fastq_file_short_r2)
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
        parser = FastqFileParser(fastq_r1=fastq_file_r1)
        line = next(parser._file_parser_handle(parser.fastq_r1))
        assert line.startswith("@") and len(line.split(":")) == 10


class TestSampleSheetParser:

    def test_get_barcode_to_sample_mapping(self, samplesheet_entries):
        parser = SampleSheetParser("this-is-a-samplesheet")
        with mock.patch.object(parser, "get_file_handle", side_effect=[samplesheet_entries]):
            mapping = list()
            for barcode, sample_id in parser.get_barcode_to_sample_mapping().items():
                mapping.append("\t".join([sample_id] + barcode.split("+")))
            assert sorted(mapping) == sorted(samplesheet_entries)

    def test_samplesheet_parser(self, samplesheet_file):
        parser = SampleSheetParser(samplesheet_file)
        mapping = parser.get_barcode_to_sample_mapping()
        expected_mapping = {
            "GGGGGGGG+AGATCTCG": "Sample1",
            "GAAGATTT+TTTACTCT": "Sample2",
            "GAAGATTT+AAAACGCC": "Sample3"}
        assert mapping == expected_mapping
