
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


class TestSampleSheetParser:

    def test_get_barcode_to_sample_mapping(self, samplesheet_entries):
        parser = SampleSheetParser("this-is-a-samplesheet")
        with mock.patch.object(parser, "get_file_handle", side_effect=[samplesheet_entries]):
            mapping = list()
            for barcode, sample_id in parser.get_barcode_to_sample_mapping().items():
                mapping.append("\t".join([sample_id] + barcode.split("+")))
            assert sorted(mapping) == sorted(samplesheet_entries)
