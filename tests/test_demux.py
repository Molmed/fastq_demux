
from unittest import mock
import os.path
import tempfile

from .context import fastq_demux
from fastq_demux.demux import FastqDemultiplexer, Demultiplexer, DemultiplexResults
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
                unknown_barcode="this-is-the-unknown-barcode-id",
                no_gzip_compression=False)
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
            unknown_barcode="Unknown",
            no_gzip_compression=True)

        expected_filename = os.path.join(
            outdir,
            f"{prefix}{sample_id}_{barcode.replace('+', '-')}_R1.fastq"
        )

        assert demux.fastq_file_name_from_sample_id(sample_id, barcode) == [expected_filename]

        fastq_parser.is_single_end = False
        demux = Demultiplexer(
            fastq_file_parser=fastq_parser,
            samplesheet_parser=samplesheet_parser,
            prefix=prefix,
            outdir=outdir,
            unknown_barcode="Unknown",
            no_gzip_compression=False)
        expected_filename = f"{expected_filename}.gz"

        assert demux.fastq_file_name_from_sample_id(sample_id, barcode) == [
            expected_filename,
            expected_filename.replace("_R1.", "_R2.")]


class TestFastqDemultiplexer:

    def test_demultiplex(self, fastq_parser, fastq_writer, fastq_records):
        demuxer = FastqDemultiplexer(fastq_parser, fastq_writer)
        results = demuxer.demultiplex()

        # ensure that the results are as expected and that the correct writers have been used
        for i, record in enumerate(fastq_records):
            expected_barcode = record[0][0].split(":")[-1]
            counter = results.known_barcodes if i > 0 else results.unknown_barcodes
            assert counter[expected_barcode] == 1
            expected_barcode = expected_barcode if i > 0 else "Unknown"
            for j, read_record in enumerate(record):
                fastq_writer[expected_barcode][j].write_record.assert_called_once_with(read_record)

    def test_single_barcode_from_record(self):
        barcode = "ACGTGT"
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


class TestDemultiplexResults:

    def test_add_functions(self):
        barcodes = [f"this-is-barcode-{i+1}" for i in range(3)]
        results = DemultiplexResults()

        counters = [results.known_barcodes, results.unknown_barcodes]
        count_fn = [results.add_known, results.add_unknown]
        for i, fn in enumerate(count_fn):
            fn(barcodes[0])
            fn(barcodes[0])
            fn(barcodes[1], n=2)
            assert counters[i][barcodes[0]] == 2
            assert counters[i][barcodes[1]] == 2
            assert counters[i][barcodes[2]] == 0

    def test_format_counts(self):
        results = DemultiplexResults()

        # add some counts to barcodes
        barcode_counts = dict(zip(
            [f"barcode-{i+1}" for i in range(10)], [7*i for i in range(10)]))
        total = sum(barcode_counts.values())
        for i, (barcode, count) in enumerate(barcode_counts.items()):
            if i < 5:
                results.add_known(barcode, n=count)
            else:
                results.add_unknown(barcode, n=count)

        # verify the summary of the known barcodes and top 2 unknown barcodes
        summarized_results = results.summarize_counts(counts=results.known_barcodes) + \
                            results.summarize_counts(counts=results.unknown_barcodes, n_values=2)
        for summary in summarized_results:
            assert summary[0] in barcode_counts
            assert summary[1] == barcode_counts[summary[0]]
            assert summary[2] == round(100. * barcode_counts[summary[0]] / total, 1)
