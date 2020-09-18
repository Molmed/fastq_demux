
import os
import pytest
import tempfile

from unittest import mock

from .context import fastq_demux
from fastq_demux.parser import FastqFileParser, SampleSheetParser
from fastq_demux.writer import FastqFileWriter, FastqFileWriterHandler
from fastq_demux.demux import DemultiplexResults


@pytest.fixture
def samplesheet_file() -> str:
    return "tests/samplesheet.tsv"


@pytest.fixture
def fastq_file_r1() -> str:
    return "tests/dual-index_Undetermined_S0_L001_R1_001.fastq.gz"


@pytest.fixture
def fastq_file_r2() -> str:
    return "tests/dual-index_Undetermined_S0_L001_R1_001.fastq.gz"


@pytest.fixture
def fastq_file_short_r2() -> str:
    return "tests/dual-index-short_Undetermined_S0_L001_R1_001.fastq.gz"


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
def fastq_writer(fastq_records):
    fastq_writer = FastqFileWriterHandler(
        prefix="this-is-a-prefix",
        outdir="this-is-the-outdir",
        is_single_end=(len(fastq_records[0]) == 1),
        no_gzip_compression=False
    )
    fastq_writer.fastq_file_writers = dict()
    for i, record in enumerate(fastq_records):
        barcode = f"record{i + 1}" if i > 0 else "Unknown"
        fastq_writer.fastq_file_writers[barcode] = [
            mock.MagicMock(spec=FastqFileWriter)() for _ in record]

    return fastq_writer


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


@pytest.fixture
def samplesheet_parser(samplesheet_entries):
    parser = SampleSheetParser(samplesheet_file="this-is-a-samplesheet-file")
    with mock.patch.object(parser, "get_file_handle", return_value=samplesheet_entries):
        yield parser


@pytest.fixture
def barcode_counts(samplesheet_parser):
    barcodes = list(samplesheet_parser.get_barcode_to_sample_mapping().keys())
    return dict(zip(barcodes, [7 * i for i, _ in enumerate(barcodes)]))


@pytest.fixture
def demultiplex_results(barcode_counts):
    results = DemultiplexResults()
    for i, (barcode, count) in enumerate(barcode_counts.items()):
        if i < 2:
            results.add_known(barcode, n=count)
        else:
            results.add_unknown(barcode, n=count)
    return results
