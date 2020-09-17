
from .context import fastq_demux
from fastq_demux.demux import FastqDemultiplexer, DemultiplexResults


class TestFastqDemultiplexer:

    def test_demultiplex(self, fastq_parser, fastq_writer, fastq_records):
        unknown_barcode = "Unknown"
        demuxer = FastqDemultiplexer(fastq_parser, fastq_writer, unknown_barcode)
        results = demuxer.demultiplex()

        # ensure that the results are as expected and that the correct writers have been used
        for i, record in enumerate(fastq_records):
            expected_barcode = record[0][0].split(":")[-1]
            counter = results.known_barcodes if i > 0 else results.unknown_barcodes
            assert counter[expected_barcode] == 1
            expected_barcode = expected_barcode if i > 0 else "Unknown"
            for j, read_record in enumerate(record):
                fastq_writer.fastq_file_writers[
                    expected_barcode][j].write_record.assert_called_once_with(read_record)

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

    def test_summarize_counts(self, barcode_counts, demultiplex_results):
        total = sum(barcode_counts.values())

        # verify the summary of the known barcodes and top 2 unknown barcodes
        summarized_results = demultiplex_results.summarize_counts(
            counts=demultiplex_results.known_barcodes) + demultiplex_results.summarize_counts(
            counts=demultiplex_results.unknown_barcodes, n_values=2)
        for summary in summarized_results:
            assert summary[0] in barcode_counts
            assert summary[1] == barcode_counts[summary[0]]
            assert summary[2] == round(100. * barcode_counts[summary[0]] / total, 1)

    def test_stats_json(self, samplesheet_parser, demultiplex_results):
        barcode_to_sample_mapping = samplesheet_parser.get_barcode_to_sample_mapping()
        stats_json = demultiplex_results.stats_json(barcode_to_sample_mapping)
        assert stats_json["UnknownBarcodes"]["Barcodes"] == demultiplex_results.unknown_barcodes
        assert stats_json["ConversionResults"]["Undetermined"]["NumberReads"] == sum(
            demultiplex_results.unknown_barcodes.values())
        for barcode_result in stats_json["ConversionResults"]["DemuxResults"]:
            barcode = barcode_result["IndexMetrics"][0]["IndexSequence"]
            assert barcode in barcode_to_sample_mapping
            assert barcode_result["IndexMetrics"][0]["MismatchCounts"]["0"] == \
                   demultiplex_results.known_barcodes[barcode]
            assert barcode_result["NumberReads"] == demultiplex_results.known_barcodes[barcode]
            assert barcode_result["SampleId"] == barcode_to_sample_mapping[barcode]
