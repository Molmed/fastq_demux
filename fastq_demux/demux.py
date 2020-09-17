import os.path

from collections import Counter
from typing import Dict, List, Optional, Tuple
from fastq_demux.parser import FastqFileParser, SampleSheetParser
from fastq_demux.writer import FastqFileWriter


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


class Demultiplexer:

    def __init__(
            self,
            fastq_file_parser: FastqFileParser,
            samplesheet_parser: SampleSheetParser,
            prefix: str,
            outdir: str,
            unknown_barcode: str,
            no_gzip_compression: bool):
        self.fastq_file_parser: FastqFileParser = fastq_file_parser
        self.samplesheet_parser: SampleSheetParser = samplesheet_parser
        self.prefix: str = prefix
        self.outdir: str = outdir
        self.unknown_barcode: str = unknown_barcode
        self.no_gzip_compression: bool = no_gzip_compression
        self.is_single_end: bool = fastq_file_parser.is_single_end
        self.barcode_fastq_writers: Dict[str, List[FastqFileWriter]] = dict()

    def get_barcode_file_writers(self) -> Dict[str, List[FastqFileWriter]]:
        barcode_to_sample_mapping: Dict[str, str] = \
            self.samplesheet_parser.get_barcode_to_sample_mapping()
        self.barcode_fastq_writers: Dict[str, List[FastqFileWriter]] = dict()
        # append an entry for the unknown barcodes
        barcode_to_sample_mapping[self.unknown_barcode] = "Sample"
        for barcode, sample_id in barcode_to_sample_mapping.items():
            self.barcode_fastq_writers[barcode] = [
                FastqFileWriter(fastq_file)
                for fastq_file in self.fastq_file_name_from_sample_id(sample_id, barcode)]
        return self.barcode_fastq_writers

    def fastq_file_name_from_sample_id(self, sample_id: str, barcode: str) -> List[str]:
        ext = "fastq.gz" if not self.no_gzip_compression else "fastq"
        return [
            os.path.join(
                self.outdir,
                f"{self.prefix}{sample_id}_{barcode.replace('+', '-')}_R{read_no}.{ext}"
            ) for read_no in range(1, 3 - int(self.is_single_end))]

    def demultiplex(self) -> DemultiplexResults:
        self.get_barcode_file_writers()
        demuxer: FastqDemultiplexer = FastqDemultiplexer(
            fastq_parser=self.fastq_file_parser,
            fastq_writer=self.barcode_fastq_writers,
            unknown_barcode=self.unknown_barcode)
        return demuxer.demultiplex()


class FastqDemultiplexer:

    def __init__(
            self,
            fastq_parser: FastqFileParser,
            fastq_writer: Dict[str, List[FastqFileWriter]],
            unknown_barcode: str = "Unknown"):
        self.fastq_parser: FastqFileParser = fastq_parser
        self.fastq_writer: Dict[str, List[FastqFileWriter]] = fastq_writer
        self.demultiplex_results: DemultiplexResults = DemultiplexResults()
        self.unknown_barcode: str = unknown_barcode

    def demultiplex(self) -> DemultiplexResults:
        for record in self.fastq_parser.fastq_records():
            self.demultiplex_record(record)
        return self.demultiplex_results

    def demultiplex_record(self, fastq_record: List[List[str]]):
        barcode: str = self.barcode_from_record(fastq_record[0])

        def _writers(bcode, count_fn) -> List[FastqFileWriter]:
            writer_list: List[FastqFileWriter] = self.fastq_writer[bcode]
            count_fn(barcode)
            return writer_list

        try:
            writers: List[FastqFileWriter] = _writers(
                barcode,
                self.demultiplex_results.add_known)
        except KeyError:
            writers: List[FastqFileWriter] = _writers(
                self.unknown_barcode,
                self.demultiplex_results.add_unknown)

        for read_no, record in enumerate(fastq_record):
            writers[read_no].write_record(record)

    @classmethod
    def barcode_from_record(cls, record: List[str]) -> str:
        return record[0].split(":")[-1]
