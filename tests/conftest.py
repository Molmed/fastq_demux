
import os
import pytest
import tempfile

from io import StringIO
from unittest import mock

from .context import fastq_demux
from fastq_demux.parser import FastqFileParserHeaderIndex, FastqFileParserI1, FastqFileParserI2, \
    SampleSheetParser
from fastq_demux.writer import FastqFileWriter, FastqFileWriterHandler
from fastq_demux.demux import DemultiplexResults


@pytest.fixture
def single_index_samplesheet_file() -> str:
    return "tests/samplesheet_single_index.tsv"


@pytest.fixture
def dual_index_samplesheet_file() -> str:
    return "tests/samplesheet_dual_index.tsv"


@pytest.fixture
def fastq_file_r1() -> str:
    return "tests/dual-index_Undetermined_S0_L001_R1_001.fastq.gz"


@pytest.fixture
def single_index_fastq_file_r1() -> str:
    return "tests/single-index_Undetermined_S0_L001_R1_001.fastq.gz"


@pytest.fixture
def fastq_file_r2() -> str:
    return "tests/dual-index_Undetermined_S0_L001_R2_001.fastq.gz"


@pytest.fixture
def fastq_file_i1() -> str:
    return "tests/dual-index_Undetermined_S0_L001_I1_001.fastq.gz"


@pytest.fixture
def fastq_file_i2() -> str:
    return "tests/dual-index_Undetermined_S0_L001_I2_001.fastq.gz"


@pytest.fixture
def fastq_file_short_r2() -> str:
    return "tests/dual-index-short_Undetermined_S0_L001_R1_001.fastq.gz"


def _mock_fastq_records(paired_end=True, dual_index=False):
    records = []
    for recordno in range(1, 11):
        barcode = f"record{recordno}i1"
        barcode = f"{barcode}+record{recordno}i2" if dual_index else barcode
        records.append((
            [
                [
                    f"R{read_no}, row{row}:{barcode}"
                    for row in range(1, 5)]
                for read_no in range(1, 2 + int(paired_end))],
            barcode
        ))
    return records


@pytest.fixture
def single_end_dual_index_fastq_records():
    return _mock_fastq_records(paired_end=False, dual_index=True)

@pytest.fixture
def single_end_single_index_fastq_records():
    return _mock_fastq_records(paired_end=False, dual_index=False)


@pytest.fixture
def paired_end_single_index_fastq_records():
    return _mock_fastq_records(paired_end=True, dual_index=False)


@pytest.fixture
def paired_end_dual_index_fastq_records():
    return _mock_fastq_records(paired_end=True, dual_index=True)


def _fastq_record_index_header(barcode):
    record = [
        ":".join([
            "@",
            "this",
            "-is-",
            "a",
            " _header_",
            "w1th[",
            "s]ome",
            "0dd",
            "characters",
            "and",
            barcode]),
        "this-is-the-nucleotide-sequence",
        "+",
        "this-is-the-quality-sequence"
    ]
    return record, barcode


@pytest.fixture
def fastq_record_dual_index_header():
    return _fastq_record_index_header("ACGTGTAG+TGACATGA")


@pytest.fixture
def fastq_record_single_index_header():
    return _fastq_record_index_header("ACGTGT")


def _mock_fastq_writer(mock_fastq_records):
    fastq_writer = FastqFileWriterHandler(
        prefix="this-is-a-prefix",
        outdir="this-is-the-outdir",
        is_single_end=(len(mock_fastq_records[0][0]) == 1),
        no_gzip_compression=False
    )
    fastq_writer.fastq_file_writers = dict()
    for i, (record, barcode) in enumerate(mock_fastq_records):
        barcode = barcode if i > 0 else "Unknown"
        fastq_writer.fastq_file_writers[barcode] = [
            mock.MagicMock(spec=FastqFileWriter)() for _ in record]

    return fastq_writer


@pytest.fixture
def single_end_single_index_fastq_writer(single_end_single_index_fastq_records):
    return _mock_fastq_writer(single_end_single_index_fastq_records)


@pytest.fixture
def single_end_dual_index_fastq_writer(single_end_dual_index_fastq_records):
    return _mock_fastq_writer(single_end_dual_index_fastq_records)


@pytest.fixture
def paired_end_single_index_fastq_writer(paired_end_single_index_fastq_records):
    return _mock_fastq_writer(paired_end_single_index_fastq_records)


@pytest.fixture
def paired_end_dual_index_fastq_writer(paired_end_dual_index_fastq_records):
    return _mock_fastq_writer(paired_end_dual_index_fastq_records)


def _mock_fastq_handles(mock_fastq_records):
    read_handles = []
    index_handles = []
    for records, barcode in mock_fastq_records:
        for idxno, idxseq in enumerate(barcode.split("+")):
            if not len(index_handles) > idxno:
                index_handles.append(StringIO())
            index_handles[idxno].write(
                "\n".join([records[0][0], idxseq, "+", "F"*len(idxseq), ""]))
            index_handles[idxno].flush()
        for readno, record in enumerate(records):
            if not len(read_handles) > readno:
                read_handles.append(StringIO())
            for row in record:
                read_handles[readno].write(f"{row}\n")
            read_handles[readno].flush()
    for handle in read_handles + index_handles:
        handle.seek(0)
    return read_handles, index_handles


@pytest.fixture
def single_end_single_index_fastq_handles(single_end_single_index_fastq_records):
    return _mock_fastq_handles(single_end_single_index_fastq_records)


@pytest.fixture
def single_end_dual_index_fastq_handles(single_end_dual_index_fastq_records):
    return _mock_fastq_handles(single_end_dual_index_fastq_records)


@pytest.fixture
def paired_end_single_index_fastq_handles(paired_end_single_index_fastq_records):
    return _mock_fastq_handles(paired_end_single_index_fastq_records)


@pytest.fixture
def paired_end_dual_index_fastq_handles(paired_end_dual_index_fastq_records):
    return _mock_fastq_handles(paired_end_dual_index_fastq_records)


def _mock_fastq_parser(mock_fastq_records, parser_class):
    def _fastq_records():
        yield from mock_fastq_records

    parser = mock.MagicMock(spec=parser_class)
    parser.return_value.fastq_records.side_effect = _fastq_records
    parser.return_value.barcode_from_record = parser_class.barcode_from_record
    return parser.return_value


@pytest.fixture
def single_end_single_index_header_fastq_parser(single_end_single_index_fastq_records):
    return _mock_fastq_parser(single_end_single_index_fastq_records, FastqFileParserHeaderIndex)


@pytest.fixture
def single_end_dual_index_header_fastq_parser(single_end_dual_index_fastq_records):
    return _mock_fastq_parser(single_end_dual_index_fastq_records, FastqFileParserHeaderIndex)


@pytest.fixture
def paired_end_single_index_header_fastq_parser(paired_end_single_index_fastq_records):
    return _mock_fastq_parser(paired_end_single_index_fastq_records, FastqFileParserHeaderIndex)


@pytest.fixture
def paired_end_dual_index_header_fastq_parser(paired_end_dual_index_fastq_records):
    return _mock_fastq_parser(paired_end_dual_index_fastq_records, FastqFileParserHeaderIndex)


@pytest.fixture
def single_end_single_index_fastq_parser(single_end_single_index_fastq_records):
    return _mock_fastq_parser(single_end_single_index_fastq_records, FastqFileParserI1)


@pytest.fixture
def paired_end_single_index_fastq_parser(paired_end_single_index_fastq_records):
    return _mock_fastq_parser(paired_end_single_index_fastq_records, FastqFileParserI1)


@pytest.fixture
def single_end_dual_index_fastq_parser(single_end_dual_index_fastq_records):
    return _mock_fastq_parser(single_end_dual_index_fastq_records, FastqFileParserI2)


@pytest.fixture
def paired_end_dual_index_fastq_parser(paired_end_dual_index_fastq_records):
    return _mock_fastq_parser(paired_end_dual_index_fastq_records, FastqFileParserI2)


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
def demultiplex_results(barcode_counts, samplesheet_parser):
    barcode_to_sample_mapping = samplesheet_parser.get_barcode_to_sample_mapping()
    # create a DemultiplexResults object using only a subset of the barcodes as "knowns"
    known_barcodes = list(barcode_to_sample_mapping.keys())[2:]
    results = DemultiplexResults(
        barcode_to_sample_mapping={
            barcode: barcode_to_sample_mapping[barcode] for barcode in known_barcodes})

    # add some counts to the results
    for barcode, count in barcode_counts.items():
        results.add(barcode, n=count)
    return results
