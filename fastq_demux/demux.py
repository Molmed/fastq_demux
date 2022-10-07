
from __future__ import annotations

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
        self.barcode_mismatches_count: Dict[str, Counter] = {
            barcode: Counter() for barcode in barcode_to_sample_mapping.keys()}

    def add(self, barcode: str, n: int = 1, distance: int = 0) -> None:
        """
        Increase the count of the supplied barcode by n (default = 1)

        :param barcode: the barcode to increase the count for
        :param n: the value to increase the count with (default = 1)
        :param distance: the number of mismatches encountered in the barcode (default = 0)
        """
        self.barcode_counts[barcode] += n
        try:
            self.barcode_mismatches_count[barcode][str(distance)] += n
        except KeyError:
            pass

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

    def _sample_to_barcode(self) -> Dict[str, List[str]]:
        sample_to_barcode = {
            sample: []
            for sample in set(list(self.barcode_to_sample_mapping.values()))}
        for barcode, sample in self.barcode_to_sample_mapping.items():
            sample_to_barcode[sample].append(barcode)
        return sample_to_barcode

    def barcode_mismatch_counts(self, barcode: str) -> List[int]:
        max_distance = max([0] + [
            int(distance)
            for distance in self.barcode_mismatches_count[barcode].keys()])
        return [
            self.barcode_mismatches_count[barcode][str(distance)]
            for distance in range(max_distance + 1)]

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

        for sample, barcodes in self._sample_to_barcode().items():
            demux_results.append({
                "SampleId": sample,
                "NumberReads": sum([self.barcode_counts[barcode] for barcode in barcodes]),
                "IndexMetrics": [{
                    "IndexSequence": barcode,
                    "MismatchCounts": {
                        str(distance): self.barcode_mismatches_count[barcode][str(distance)]
                        for distance, count in enumerate(self.barcode_mismatch_counts(barcode))
                    }
                } for barcode in barcodes]
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
            unknown_barcode: str,
            mismatches: int = 0):
        """
        Create a new FastqDemultiplexer instance

        :param fastq_parser: a FastqFileParser instance used to parse the input FASTQ file(s)
        :param fastq_writer: a FastqFileWriterHandler instance used to write the demultiplexed
        output to corresponding FASTQ file(s)
        :param demultiplex_results: a Demultiplexresults instance used for collecting statistics
        :param unknown_barcode: a string to collect the unknown barcodes under
        :param mismatches: the number of mismatches allowed when matching barcodes (default 0)
        """
        self.fastq_parser: FastqFileParser = fastq_parser
        self.fastq_writer: FastqFileWriterHandler = fastq_writer
        self.demultiplex_results: DemultiplexResults = demultiplex_results
        self.unknown_barcode: str = unknown_barcode
        self.mismatches = mismatches
        self.mismatched_barcodes: Dict[str, Tuple[str, int, str]] = {}

    @staticmethod
    def create_fastq_demultiplexer(
        fastq_parser: FastqFileParser,
        fastq_writer: FastqFileWriterHandler,
        demultiplex_results: DemultiplexResults,
        unknown_barcode: str,
        mismatches: int = 0) -> FastqDemultiplexer:
        if mismatches == 0:
            return FastqDemultiplexer(
                fastq_parser,
                fastq_writer,
                demultiplex_results,
                unknown_barcode,
                mismatches)
        return FastqMismatchDemultiplexer(
                fastq_parser,
                fastq_writer,
                demultiplex_results,
                unknown_barcode,
                mismatches)

    def demultiplex(self) -> DemultiplexResults:
        """
        Iterate over the fastq_parser and write output to the fastq_writer

        :return: the DemultiplexResults instance containing the statistics from the demultiplexing
        """
        for record, barcode in self.fastq_parser.fastq_records():
            self.demultiplex_record(record, barcode)
        return self.demultiplex_results

    def _write_matching_barcode(self, fastq_record: List[List[str]], barcode: str) -> None:
        self.fastq_writer.write_fastq_record(barcode, fastq_record)

    def demultiplex_record(self, fastq_record: List[List[str]], barcode: str) -> None:
        """
        Demultiplex a FASTQ read (or read pair) using the barcode in the FASTQ header and write it
        to the corresponding FastqFileWriter

        :param fastq_record: a 1- or 2-element list of FASTQ reads for single-end or paired-end,
        respectively. Each read is a list of 4 strings
        :param barcode: a string containing the barcode belonging to the record
        """
        try:
            self._write_matching_barcode(fastq_record, barcode)
        except KeyError:
            self._write_matching_barcode(fastq_record, self.unknown_barcode)
        self.demultiplex_results.add(barcode)


class FastqMismatchDemultiplexer(FastqDemultiplexer):

    @staticmethod
    def hamming_distance(str1: str, str2: str) -> int:
        return sum([int(s1 != s2) for (s1, s2) in zip(str1, str2)])

    def match_mismatched_single_barcode(
            self, barcode: str, known_barcodes: List[str]) -> Tuple[str, int]:
        distances = list(map(
            lambda x: self.hamming_distance(x, barcode),
            known_barcodes))
        no_mm = -1
        while no_mm < self.mismatches:
            no_mm += 1
            try:
                i = distances.index(no_mm)
                if distances.count(no_mm) > 1 and \
                        len(set([
                            known_barcodes[x] for x, d in enumerate(distances) if d == no_mm])) > 1:
                    raise Exception(
                        f"the barcode {barcode} is ambiguous when allowing {self.mismatches} "
                        f"mismatches")
                return known_barcodes[i], no_mm
            except ValueError:
                pass
        return "", -1

    def match_mismatched_barcode(self, barcode: str) -> None:
        keys = list(
            self.demultiplex_results.barcode_to_sample_mapping.keys())
        matched_key = []
        distance = 0
        for ix, sindex in enumerate(barcode.split("+")):
            known_barcodes = [k.split("+")[ix] for k in keys if k != self.unknown_barcode]
            mk, d = self.match_mismatched_single_barcode(sindex, known_barcodes)
            matched_key.append(mk)
            distance = max([d, distance])
        matched_barcode = "+".join(matched_key)
        counted_barcode = matched_barcode
        if matched_barcode not in keys:
            matched_barcode = self.unknown_barcode
            distance = 0
            counted_barcode = barcode
        self.mismatched_barcodes[barcode] = (matched_barcode, distance, counted_barcode)

    def demultiplex_record(self, fastq_record: List[List[str]], barcode: str) -> None:
        """
        Demultiplex a FASTQ read (or read pair) using the barcode in the FASTQ header and write it
        to the corresponding FastqFileWriter

        :param fastq_record: a 1- or 2-element list of FASTQ reads for single-end or paired-end,
        respectively. Each read is a list of 4 strings
        :param barcode: a string containing the barcode belonging to the record
        """
        try:
            self._write_matching_barcode(fastq_record, self.mismatched_barcodes[barcode][0])
        except KeyError:
            self.match_mismatched_barcode(barcode)
            return self.demultiplex_record(fastq_record, barcode)
        self.demultiplex_results.add(
            self.mismatched_barcodes[barcode][2],
            distance=self.mismatched_barcodes[barcode][1])
