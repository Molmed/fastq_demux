
import os

from typing import Dict, List, Optional, TextIO
from fastq_demux.parser import FastqFileParser


class FastqFileWriter:

    def __init__(self, fastq_file: str):
        self.fastq_file: str = fastq_file
        self.fastq_file_handle: Optional[TextIO] = None
        self.open_handle()

    def open_handle(self) -> None:
        if self.fastq_file_handle is None:
            self.fastq_file_handle = FastqFileParser.open_func(self.fastq_file)(
                self.fastq_file, "wt")

    def close_handle(self) -> None:
        if self.fastq_file_handle is not None:
            self.fastq_file_handle.close()

    def write_record(self, record: List[str]) -> None:
        self.fastq_file_handle.write("\n".join(record + [""]))


class FastqFileWriterHandler:

    def __init__(
            self,
            prefix: str,
            outdir: str,
            is_single_end: bool,
            no_gzip_compression: bool):

        self.prefix: str = prefix
        self.outdir: str = outdir
        self.is_single_end: bool = is_single_end
        self.no_gzip_compression: bool = no_gzip_compression

        self.fastq_file_writers: Optional[Dict[str, List[FastqFileWriter]]] = None

    def fastq_file_name_from_sample_id(
            self, sample_id: str, barcode: str) -> List[str]:
        ext: str = "fastq.gz" if not self.no_gzip_compression else "fastq"
        return [
            os.path.join(
                self.outdir,
                f"{self.prefix}{sample_id}_{barcode.replace('+', '-')}_R{read_no}.{ext}"
            ) for read_no in range(1, 3 - int(self.is_single_end))]

    def fastq_file_writers_from_sample_id(
            self, sample_id: str, barcode: str) -> List[FastqFileWriter]:
        return [
            FastqFileWriter(fastq_file)
            for fastq_file in self.fastq_file_name_from_sample_id(
                sample_id,
                barcode)]

    def fastq_file_writers_from_mapping(
            self,
            barcode_to_sample_mapping: Dict[str, str]) -> Dict[str, List[FastqFileWriter]]:

        if not self.fastq_file_writers:
            self.fastq_file_writers: Dict[str, List[FastqFileWriter]] = dict()
            for barcode, sample_id in barcode_to_sample_mapping.items():
                self.fastq_file_writers[barcode] = \
                    self.fastq_file_writers_from_sample_id(sample_id, barcode)

        return self.fastq_file_writers

    def write_fastq_record(self, barcode: str, fastq_record: List[List[str]]) -> None:
        for read_no, record in enumerate(fastq_record):
            self.fastq_file_writers[barcode][read_no].write_record(record)
