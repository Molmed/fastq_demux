
import gzip
from typing import Dict, Iterator, List, Optional, TextIO


class FastqFileParser:

    def __init__(self, fastq_r1: str, fastq_r2: Optional[str] = None):
        self.fastq_r1: str = fastq_r1
        self.fastq_r2: Optional[str] = fastq_r2
        self.is_single_end: bool = self.fastq_r2 is None

    @classmethod
    def _ensure_generator_empty(cls, gen: Iterator) -> bool:
        try:
            next(gen)
        except StopIteration:
            return True
        return False

    def ensure_filehandles_empty(self, handles: List[TextIO]) -> bool:
        """
        Checks whether supplied filehandles are empty. Note that this will consume entries from the handles.
        """
        return all([self._ensure_generator_empty(handle) for handle in handles])

    def fastq_records(self) -> Iterator[List[List[str]]]:
        handles: List[TextIO] = self.get_file_handles()
        try:
            while True:
                yield [[next(handle).strip() for _ in range(4)] for handle in handles]
        except StopIteration:
            pass
        if not self.ensure_filehandles_empty(handles):
            raise Exception(f"Fastq files {self.fastq_r1} and {self.fastq_r2} are not of the same length!")

    def get_file_handles(self) -> List[TextIO]:
        handles: List[TextIO] = [self._gzfile_parser(self.fastq_r1)]
        if not self.is_single_end:
            handles.append(self._gzfile_parser(self.fastq_r2))
        return handles

    @classmethod
    def _gzfile_parser(cls, gzipped_file: str) -> TextIO:
        with gzip.open(gzipped_file, 'rt') as fh:
            yield from fh


class SampleSheetParser:

    def __init__(self, samplesheet_file: str):
        self.samplesheet_file: str = samplesheet_file

    def get_file_handle(self) -> TextIO:
        with open(self.samplesheet_file, 'rt') as fh:
            yield from fh

    def get_barcode_to_sample_mapping(self) -> Dict[str, str]:
        barcode_to_sample: Dict[str, str] = dict()
        for row in self.get_file_handle():
            pcs: List[str] = row.strip().split("\t")
            barcode_to_sample["+".join(pcs[1:])] = pcs[0]
        return barcode_to_sample
