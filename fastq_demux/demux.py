import os.path

from collections import Counter
from typing import Dict, List, Optional, Tuple
from fastq_demux.parser import FastqFileParser, SampleSheetParser
from fastq_demux.writer import FastqFileWriter


class Demultiplexer:

    def __init__(
            self,
            fastq_file_r1: str,
            samplesheet: str,
            prefix: str,
            outdir: str,
            unknown_barcode: str,
            fastq_file_r2: Optional[str] = None):
        self.fastq_file_r1: str = fastq_file_r1
        self.fastq_file_r2: Optional[str] = fastq_file_r2
        self.samplesheet: str = samplesheet
        self.prefix: str = prefix
        self.outdir: str = outdir
        self.unknown_barcode = unknown_barcode
        self.is_single_end = fastq_file_r2 is None

    def get_barcode_file_writers(self) -> Dict[str, List[FastqFileWriter]]:
        parser: SampleSheetParser = SampleSheetParser(self.samplesheet)
        barcode_to_sample_mapping: Dict[str, str] = parser.get_barcode_to_sample_mapping()
        barcode_to_file_writer: Dict[str, List[FastqFileWriter]] = dict()
        # append an entry for the unknown barcodes
        barcode_to_sample_mapping[self.unknown_barcode] = "Sample"
        for barcode, sample_id in barcode_to_sample_mapping.items():
            barcode_to_file_writer[barcode] = [
                FastqFileWriter(fastq_file) for fastq_file in self.fastq_file_name_from_sample_id(sample_id, barcode)]
        return barcode_to_file_writer

    def fastq_file_name_from_sample_id(self, sample_id: str, barcode: str) -> List[str]:
        return [
            os.path.join(
                self.outdir,
                f"{self.prefix}{sample_id}_{barcode.replace('+', '-')}_R{read_no}.fastq.gz"
            ) for read_no in range(1, 3 - int(self.is_single_end))]

    def demultiplex(self) -> Tuple[Counter, Counter]:
        barcode_writers: Dict[str, List[FastqFileWriter]] = self.get_barcode_file_writers()
        parser: FastqFileParser = FastqFileParser(self.fastq_file_r1, self.fastq_file_r2)
        demuxer: FastqDemultiplexer = FastqDemultiplexer(
            fastq_parser=parser,
            fastq_writer=barcode_writers,
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
