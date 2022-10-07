
import pytest

from .context import fastq_demux
from fastq_demux.demux import FastqDemultiplexer, FastqMismatchDemultiplexer, DemultiplexResults


class TestFastqDemultiplexer:

    @staticmethod
    def _helper_demultiplex(parser, writer, records):
        unknown_barcode = "Unknown"
        demultiplex_results = DemultiplexResults(barcode_to_sample_mapping={})
        demuxer = FastqDemultiplexer.create_fastq_demultiplexer(
            parser, writer, demultiplex_results, unknown_barcode, mismatches=0)
        demuxer.demultiplex()

        # ensure that the results are as expected and that the correct writers have been used
        for i, (record, expected_barcode) in enumerate(records):
            assert demultiplex_results.barcode_counts[expected_barcode] == 1
            expected_barcode = expected_barcode if i > 0 else "Unknown"
            for j, read_record in enumerate(record):
                writer.fastq_file_writers[
                    expected_barcode][j].write_record.assert_called_once_with(read_record)

    def test_demultiplex_single_end_index(
            self,
            single_end_dual_index_fastq_parser,
            single_end_dual_index_fastq_writer,
            single_end_dual_index_fastq_records):
        self._helper_demultiplex(
            single_end_dual_index_fastq_parser,
            single_end_dual_index_fastq_writer,
            single_end_dual_index_fastq_records)

    def test_demultiplex_header_index(
            self,
            single_end_dual_index_header_fastq_parser,
            single_end_dual_index_fastq_writer,
            single_end_dual_index_fastq_records):
        self._helper_demultiplex(
            single_end_dual_index_header_fastq_parser,
            single_end_dual_index_fastq_writer,
            single_end_dual_index_fastq_records)

    def test_demultiplex_single_index(
            self,
            single_end_single_index_fastq_parser,
            single_end_single_index_fastq_writer,
            single_end_single_index_fastq_records):
        self._helper_demultiplex(
            single_end_single_index_fastq_parser,
            single_end_single_index_fastq_writer,
            single_end_single_index_fastq_records)

    def test_demultiplex_dual_index(
            self,
            paired_end_dual_index_fastq_parser,
            paired_end_dual_index_fastq_writer,
            paired_end_dual_index_fastq_records):
        self._helper_demultiplex(
            paired_end_dual_index_fastq_parser,
            paired_end_dual_index_fastq_writer,
            paired_end_dual_index_fastq_records)


class TestFastqMismatchDemultiplexer:

    @staticmethod
    def _get_demultiplex_results(known_barcodes):
        return DemultiplexResults(
            barcode_to_sample_mapping={
                bc: f"sample_{bc}" for bc in known_barcodes
            })

    def test_hamming_distance(self):
        assert FastqMismatchDemultiplexer.hamming_distance("ABCD", "ABCD") == 0
        assert FastqMismatchDemultiplexer.hamming_distance("ABCD", "ABC_") == 1
        assert FastqMismatchDemultiplexer.hamming_distance("ABCD", "__CD") == 2
        # it is expected that it's only along the length of the shortest string that the comparison
        # is made
        assert FastqMismatchDemultiplexer.hamming_distance("ABCD", "ABC") == 0

    def test_match_mismatched_single_barcode(
            self,
            fastq_demuxer):
        known_barcodes = [
            "AAAAAA",
            "BBBBBB",
            "CCCCCC"
            "BCBCBC"
        ]
        assert fastq_demuxer.match_mismatched_single_barcode(
            "BBBBBB",
            known_barcodes) == ("BBBBBB", 0)
        assert fastq_demuxer.match_mismatched_single_barcode(
            "BBBCBB",
            known_barcodes) == ("BBBBBB", 1)
        assert fastq_demuxer.match_mismatched_single_barcode(
            "BBCCBB",
            known_barcodes) == ("", -1)
        fastq_demuxer.mismatches = 2
        assert fastq_demuxer.match_mismatched_single_barcode(
            "BBCCBB",
            known_barcodes) == ("BBBBBB", 2)
        with pytest.raises(Exception):
            fastq_demuxer.mismatches = 3
            fastq_demuxer.match_mismatched_single_barcode(
                "CBCBCB",
                known_barcodes)

    def test_single_index_match_mismatched_barcode(
            self,
            fastq_demuxer):
        known_barcodes = [
            "AAAAAA",
            "BBBBBB",
            "CCCCCC"
            "BCBCBC"
        ]
        expected_mismatched_barcodes = {
            "BBBBBB": ("BBBBBB", 0, "BBBBBB"),
            "BBBCBB": ("BBBBBB", 1, "BBBBBB"),
            "BBCCBB": (fastq_demuxer.unknown_barcode, 0, "BBCCBB")
        }
        fastq_demuxer.demultiplex_results = self._get_demultiplex_results(known_barcodes)
        for barcode in expected_mismatched_barcodes.keys():
            fastq_demuxer.match_mismatched_barcode(barcode)
        assert fastq_demuxer.mismatched_barcodes == expected_mismatched_barcodes

    def test_dual_index_match_mismatched_barcode(
            self,
            fastq_demuxer):
        known_barcodes = [
            "AAAAAA+BBBBBB",
            "BBBBBB+CCCCCC",
            "BCBCBC+CBCBCB"
        ]
        expected_mismatched_barcodes = {
            "BBBBBB+CCCCCC": ("BBBBBB+CCCCCC", 0, "BBBBBB+CCCCCC"),
            "BBBCBB+CCCCCC": ("BBBBBB+CCCCCC", 1, "BBBBBB+CCCCCC"),
            "BBCCBB+CCCCCC": (fastq_demuxer.unknown_barcode, 0, "BBCCBB+CCCCCC"),
            "BCBCBC+CBCCCB": ("BCBCBC+CBCBCB", 1, "BCBCBC+CBCBCB"),
            "BCBBBC+CBCCCB": ("BCBCBC+CBCBCB", 1, "BCBCBC+CBCBCB"),
            "CBCBCB+BCBCBC": (fastq_demuxer.unknown_barcode, 0, "CBCBCB+BCBCBC")
        }
        fastq_demuxer.demultiplex_results = self._get_demultiplex_results(known_barcodes)
        for barcode in expected_mismatched_barcodes.keys():
            fastq_demuxer.match_mismatched_barcode(barcode)
        assert fastq_demuxer.mismatched_barcodes == expected_mismatched_barcodes

    def test_demultiplex_record(self, fastq_demuxer):
        fastq_record = [[""]]
        known_barcodes = [
            "AAAAAA",
            "BBBBBB",
            "CCCCCC",
            "BCBCBC"
        ]
        expected_mismatched_barcodes = {
            "BBBBBB": ("BBBBBB", 0, "BBBBBB"),
            "BBBCBB": ("BBBBBB", 1, "BBBBBB"),
            "BBCCBB": (fastq_demuxer.unknown_barcode, 0, "BBCCBB")
        }
        fastq_demuxer.demultiplex_results = self._get_demultiplex_results(known_barcodes)
        for test_barcode, expected_barcode in zip(
                ["BBBBBB", "BBBCBB", "BBCCBB"],
                [known_barcodes[1], known_barcodes[1], fastq_demuxer.unknown_barcode]):
            fastq_demuxer.demultiplex_record(fastq_record, test_barcode)
            fastq_demuxer._write_matching_barcode.assert_called_once_with(
                fastq_record,
                expected_barcode)
            fastq_demuxer._write_matching_barcode.reset_mock()
        assert fastq_demuxer.mismatched_barcodes == expected_mismatched_barcodes
        assert fastq_demuxer.demultiplex_results.summarize_counts() == [
            ("BBBBBB", 2, 66.7),
            ("BBCCBB", 1, 33.3)]


class TestDemultiplexResults:

    def test_add(self):
        barcodes = [f"this-is-barcode-{i+1}" for i in range(2)]
        barcode_to_sample_mapping = dict(zip(barcodes, ["" for _ in barcodes]))
        results = DemultiplexResults(barcode_to_sample_mapping=barcode_to_sample_mapping)

        results.add(barcodes[0])
        results.add(barcodes[0])
        results.add(barcodes[1], n=5)
        assert results.barcode_counts[barcodes[0]] == 2
        assert results.barcode_counts[barcodes[1]] == 5

    def test_summarize_counts(self, barcode_counts, demultiplex_results):
        total = sum(barcode_counts.values())

        # verify the summary of the known barcodes and top 2 unknown barcodes
        summarized_results = demultiplex_results.summarize_counts(barcode_set="known") + \
                             demultiplex_results.summarize_counts(barcode_set="unknown", n_values=2)
        for summary in summarized_results:
            assert summary[0] in barcode_counts
            assert summary[1] == barcode_counts[summary[0]]
            assert summary[2] == round(100. * barcode_counts[summary[0]] / total, 1)

    def test_stats_json(self, samplesheet_parser, demultiplex_results):
        barcode_to_sample_mapping = samplesheet_parser.get_barcode_to_sample_mapping()
        stats_json = demultiplex_results.stats_json()
        assert stats_json["UnknownBarcodes"]["Barcodes"] == \
               demultiplex_results.unknown_barcode_counts()
        assert stats_json["ConversionResults"]["Undetermined"]["NumberReads"] == sum(
            demultiplex_results.unknown_barcode_counts().values())
        for barcode_result in stats_json["ConversionResults"]["DemuxResults"]:
            barcode = barcode_result["IndexMetrics"][0]["IndexSequence"]
            assert barcode in barcode_to_sample_mapping
            assert barcode_result["IndexMetrics"][0]["MismatchCounts"]["0"] == \
                   demultiplex_results.known_barcode_counts()[barcode]
            assert barcode_result["NumberReads"] == \
                   demultiplex_results.known_barcode_counts()[barcode]
            assert barcode_result["SampleId"] == barcode_to_sample_mapping[barcode]
