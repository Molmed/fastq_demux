
import pytest

from io import StringIO
from .context import fastq_demux
from fastq_demux.writer import FastqFileWriter
from fastq_demux.parser import FastqFileParser


class TestFastqFileWriter:

    def test_write_record(self, fastq_records, fastq_output_file_r1):
        writer = FastqFileWriter(fastq_output_file_r1)
        for record in fastq_records:
            writer.write_record(record[0])
        writer.close_handle()
        parser = FastqFileParser(fastq_output_file_r1)
        assert list(parser.fastq_records()) == [[record[0]] for record in fastq_records]
