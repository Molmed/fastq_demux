import os.path

from collections import Counter
from typing import Dict, List, Optional, Tuple
from fastq_demux.parser import FastqFileParser, SampleSheetParser
from fastq_demux.writer import FastqFileWriter


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

    def demultiplex(self) -> Tuple[Counter, Counter]:
        self.get_barcode_file_writers()
        demuxer: FastqDemultiplexer = FastqDemultiplexer(
            fastq_parser=self.fastq_file_parser,
            fastq_writer=self.barcode_fastq_writers,
            unknown_barcode=self.unknown_barcode)
        return demuxer.demultiplex()

    @classmethod
    def format_counts(
            cls, counts: Counter, n_values: Optional[int] = None) -> List[Tuple[str, int, float]]:
        return [(
            barcode,
            count,
            round(100. * count / sum(counts.values()), 1))
            for barcode, count in counts.most_common(n_values)]


class FastqDemultiplexer:

    def __init__(
            self,
            fastq_parser: FastqFileParser,
            fastq_writer: Dict[str, List[FastqFileWriter]],
            unknown_barcode: str = "Unknown"):
        self.fastq_parser: FastqFileParser = fastq_parser
        self.fastq_writer: Dict[str, List[FastqFileWriter]] = fastq_writer
        self.unknown_barcode: str = unknown_barcode
        self.known_barcode_counts: Counter = Counter()
        self.unknown_barcode_counts: Counter = Counter()

    def demultiplex(self) -> Tuple[Counter, Counter]:
        for record in self.fastq_parser.fastq_records():
            self.demultiplex_record(record)
        return self.known_barcode_counts, self.unknown_barcode_counts

    def demultiplex_record(self, fastq_record: List[List[str]]):
        barcode: str = self.barcode_from_record(fastq_record[0])
        try:
            writers: List[FastqFileWriter] = self.fastq_writer[barcode]
            self.known_barcode_counts[barcode] += 1
        except KeyError:
            writers: List[FastqFileWriter] = self.fastq_writer[self.unknown_barcode]
            self.unknown_barcode_counts[barcode] += 1

        for read_no, record in enumerate(fastq_record):
            writers[read_no].write_record(record)

    @classmethod
    def barcode_from_record(cls, record: List[str]) -> str:
        return record[0].split(":")[-1]
