
import gzip

from typing import List, Optional, TextIO


class FastqFileWriter:

    def __init__(self, fastq_file: str):
        self.fastq_file: str = fastq_file
        self.fastq_file_handle: Optional[TextIO] = None
        self.open_handle()

    def open_handle(self) -> None:
        if self.fastq_file_handle is None:
            self.fastq_file_handle = gzip.open(self.fastq_file, "wt")

    def close_handle(self) -> None:
        if self.fastq_file_handle is not None:
            self.fastq_file_handle.close()

    def write_record(self, record: List[str]) -> None:
        self.fastq_file_handle.write("\n".join(record + [""]))
