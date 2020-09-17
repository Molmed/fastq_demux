
from collections import Counter
from typing import Dict, List, Optional, Tuple
from fastq_demux.parser import FastqFileParser
from fastq_demux.writer import FastqFileWriterHandler


class DemultiplexResults:

    def __init__(self):
        self.known_barcodes: Counter = Counter()
        self.unknown_barcodes: Counter = Counter()

    def add_known(self, barcode: str, n: int = 1) -> None:
        self.known_barcodes[barcode] += n

    def add_unknown(self, barcode: str, n: int = 1) -> None:
        self.unknown_barcodes[barcode] += n

    def summarize_counts(
            self, counts: Counter, n_values: Optional[int] = None) -> List[Tuple[str, int, float]]:
        total = sum(self.known_barcodes.values()) + sum(self.unknown_barcodes.values())
        return [(
            barcode,
            count,
            round(100. * count / total, 1))
            for barcode, count in counts.most_common(n_values)]

    def stats_json(self, barcode_to_sample_mapping: Dict[str, str]) -> Dict[str, dict]:
        demux_results: List[dict] = []
        for barcode, count in self.known_barcodes.items():
            demux_results.append({
                "SampleId": barcode_to_sample_mapping[barcode],
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
                    "NumberReads": sum(self.unknown_barcodes.values())}},
            "UnknownBarcodes": {
                "Barcodes": self.unknown_barcodes}}


class FastqDemultiplexer:

    def __init__(
            self,
            fastq_parser: FastqFileParser,
            fastq_writer: FastqFileWriterHandler,
            unknown_barcode: str):
        self.fastq_parser: FastqFileParser = fastq_parser
        self.fastq_writer: FastqFileWriterHandler = fastq_writer
        self.unknown_barcode: str = unknown_barcode
        self.demultiplex_results: DemultiplexResults = DemultiplexResults()

    def demultiplex(self) -> DemultiplexResults:
        for record in self.fastq_parser.fastq_records():
            self.demultiplex_record(record)
        return self.demultiplex_results

    def demultiplex_record(self, fastq_record: List[List[str]]):
        barcode: str = self.barcode_from_record(fastq_record[0])
        try:
            self.fastq_writer.write_fastq_record(barcode, fastq_record)
            self.demultiplex_results.add_known(barcode)
        except KeyError:
            self.fastq_writer.write_fastq_record(self.unknown_barcode, fastq_record)
            self.demultiplex_results.add_unknown(barcode)

    @classmethod
    def barcode_from_record(cls, record: List[str]) -> str:
        return record[0].split(":")[-1]
