
import os
import pytest
import tempfile

from unittest import mock

from .context import fastq_demux
from fastq_demux.parser import FastqFileParser


@pytest.fixture
def fastq_file_r1() -> str:
    return "./dual-index_Undetermined_S0_L001_R1_001.fastq.gz"


@pytest.fixture
def fastq_file_r2() -> str:
    return "./dual-index_Undetermined_S0_L001_R1_001.fastq.gz"


@pytest.fixture
def fastq_file_short_r2() -> str:
    return "./dual-index-short_Undetermined_S0_L001_R1_001.fastq.gz"


@pytest.fixture
def fastq_records():
    return [
        [
            [
                f"R{read_no}, row{row}:record{record}" for row in range(1, 5)
            ] for read_no in range(1, 3)
        ] for record in range(1, 11)
    ]


@pytest.fixture
def fastq_writer():

    def _write_record(barcode, records):
        assert barcode in [record[1][0].split(":")[-1] for record in fastq_records]
        assert records in fastq_records

    writer = mock.MagicMock()
    return writer.return_value


@pytest.fixture
def fastq_parser(fastq_records):

    def _fastq_records():
        for record in fastq_records:
            yield record

    parser = mock.MagicMock(spec=FastqFileParser)
    parser.return_value.fastq_records.side_effect = _fastq_records
    return parser.return_value


@pytest.fixture
def fastq_output_file_r1():
    with tempfile.TemporaryDirectory(prefix="TestFastqFileWriter_") as tmpdir:
        fastq_file = os.path.join(tmpdir, "test_write_record.fastq.gz")
        yield fastq_file


@pytest.fixture
def samplesheet_entries():
    return [
        "\t".join(["Sample_1", "AAAAAACC", "TTTTTTGG"]),
        "\t".join(["Sample_2", "AAAACCCC", "TTTTGGGG"]),
        "\t".join(["Sample_3", "AACCCCCC", "TTGGGGGG"]),
        "\t".join(["Sample_4", "CCCCGGGG"])]
