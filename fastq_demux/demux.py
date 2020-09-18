
from collections import Counter
from typing import Dict, List, Optional, Tuple
from fastq_demux.parser import FastqFileParser
from fastq_demux.writer import FastqFileWriterHandler


class DemultiplexResults:
    """
    Class to handle the statistics from the demultiplexing
    """
    def __init__(self, barcode_to_sample_mapping: Dict[str, str]):
        """
        Create a new DemultiplexResults instance

        :param barcode_to_sample_mapping: a dict having barcodes as keys and the corresponding
        sample-id as value, representing the known barcodes used for demultiplexing
        """
        self.barcode_to_sample_mapping: Dict[str, str] = barcode_to_sample_mapping
        self.barcode_counts: Counter = Counter()

    def add(self, barcode: str, n: int = 1) -> None:
        """
        Increase the count of the supplied barcode by n (default = 1)

        :param barcode: the barcode to increase the count for
        :param n: the value to increase the count with (default = 1)
        """
        self.barcode_counts[barcode] += n

    def known_barcode_counts(self) -> Counter:
        """
        Get a Counter representing only the known barcodes
        :return: a Counter representing only the known barcodes
        """
        return Counter({barcode: count for barcode, count in self.barcode_counts.items()
                        if barcode in self.barcode_to_sample_mapping})

    def unknown_barcode_counts(self) -> Counter:
        """
        Get a Counter representing only the unknown barcodes
        :return: a Counter representing only the unknown barcodes
        """
        return Counter({barcode: count for barcode, count in self.barcode_counts.items()
                        if barcode not in self.barcode_to_sample_mapping})

    def summarize_counts(
            self,
            barcode_set: str = "all",
            n_values: Optional[int] = None) -> List[Tuple[str, int, float]]:
        """
        Summarize the counts for a set of barcodes, presenting counts and percent of total for each
        barcode in the set

        :param barcode_set: the set of barcodes to present a summary for, this must be one of
        "all", "known" or "unknown"
        :param n_values: if specified, present only the n_values top values
        :return: a list of Tuples, each element being a barcode, count and percent of total
        """
        total: int = sum(self.barcode_counts.values())
        counter: Counter = self.barcode_counts
        if barcode_set == "known":
            counter: Counter = self.known_barcode_counts()
        elif barcode_set == "unknown":
            counter: Counter = self.unknown_barcode_counts()

        return [(
            barcode,
            count,
            round(100. * count / total, 1))
            for barcode, count in counter.most_common(n_values)]

    def stats_json(self) -> Dict[str, dict]:
        """
        Summarize the results in a data structure mimicking the output in bcl2fastq's Stats.json

        :return: a dict mimicking the output in bcl2fastq's Stats.json
        """
        demux_results: List[dict] = []
        for barcode, count in self.known_barcode_counts().items():
            demux_results.append({
                "SampleId": self.barcode_to_sample_mapping[barcode],
                "NumberReads": count,
                "IndexMetrics": [
                    {
                        "IndexSequence": barcode,
                        "MismatchCounts": {
                            "0": count
                        }
                    }
                ]
            })

        return {
            "ConversionResults": {
                "DemuxResults": demux_results,
                "Undetermined": {
                    "NumberReads": sum(self.unknown_barcode_counts().values())}},
            "UnknownBarcodes": {
                "Barcodes": self.unknown_barcode_counts()}}


class FastqDemultiplexer:
    """
    Class handling the demultiplexing of a FASTQ read according to the barcode in the header
    """

    def __init__(
            self,
            fastq_parser: FastqFileParser,
            fastq_writer: FastqFileWriterHandler,
            demultiplex_results: DemultiplexResults,
            unknown_barcode: str):
        """
        Create a new FastqDemultiplexer instance

        :param fastq_parser: a FastqFileParser instance used to parse the input FASTQ file(s)
        :param fastq_writer: a FastqFileWriterHandler instance used to write the demultiplexed
        output to corresponding FASTQ file(s)
        :param demultiplex_results: a Demultiplexresults instance used for collecting statistics
        :param unknown_barcode: a string to collect the unknown barcodes under
        """
        self.fastq_parser: FastqFileParser = fastq_parser
        self.fastq_writer: FastqFileWriterHandler = fastq_writer
        self.demultiplex_results: DemultiplexResults = demultiplex_results
        self.unknown_barcode: str = unknown_barcode

    def demultiplex(self) -> DemultiplexResults:
        """
        Iterate over the fastq_parser and write output to the fastq_writer

        :return: the DemultiplexResults instance containing the statistics from the demultiplexing
        """
        for record in self.fastq_parser.fastq_records():
            self.demultiplex_record(record)
        return self.demultiplex_results

    def demultiplex_record(self, fastq_record: List[List[str]]):
        """
        Demultiplex a FASTQ read (or read pair) using the barcode in the FASTQ header and write it
        to the corresponding FastqFileWriter

        :param fastq_record: a 1- or 2-element list of FASTQ reads for single-end or paired-end,
        respectively. Each read is a list of 4 strings
        """
        barcode: str = self.barcode_from_record(fastq_record[0])
        try:
            self.fastq_writer.write_fastq_record(barcode, fastq_record)
        except KeyError:
            self.fastq_writer.write_fastq_record(self.unknown_barcode, fastq_record)
        self.demultiplex_results.add(barcode)

    @classmethod
    def barcode_from_record(cls, record: List[str]) -> str:
        """
        Extract the barcode from a FASTQ read header. The header is expected to be ":"-delimited
        according to the bcl2fastq output format and the barcode is expected to be the last element
        when splitting on ":"

        :param record: a list of strings representing a FASTQ read. Only the first element (the
        header string) will be accessed
        :return:
        """
        return record[0].split(":")[-1]
