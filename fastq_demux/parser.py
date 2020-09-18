
import gzip
from typing import Callable, Dict, Iterator, List, Optional, TextIO


class FastqFileParser:
    """
    Class that handles parsing of a single-end FASTQ file or a paired-end FASTQ file pair
    """

    def __init__(self, fastq_r1: str, fastq_r2: Optional[str] = None):
        """
        Create a new FastqFileParser instance

        :param fastq_r1: the path to the FASTQ file with forward reads (R1)
        :param fastq_r2: for paired-end, the path to the FASTQ file with reverse reads (R2)
        """
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
        Checks whether supplied filehandles are empty. Note that this will consume an entry from
        non-empty handles in the process

        :param handles: a list of file handle objects to check
        :return: True if all supplied file handles are empty, False otherwise
        """
        return all([self._ensure_generator_empty(handle) for handle in handles])

    def fastq_records(self) -> Iterator[List[List[str]]]:
        """
        Iterates over the FASTQ files

        :return: a generator which will yield FASTQ records. Each yielded element will be a list of
        length 1 for single-end and 2 for paired-end. Each list element will in turn be a list of
        4 strings, representing the 4 lines from the FASTQ file making up the FASTQ read
        :raises Exception: if, for paired-end, the supplied FASTQ files don't contain the same
        number of reads
        """
        handles: List[TextIO] = self.get_file_handles()
        try:
            while True:
                yield [[next(handle).strip() for _ in range(4)] for handle in handles]
        except StopIteration:
            pass
        if not self.ensure_filehandles_empty(handles):
            raise Exception(
                f"Fastq files {self.fastq_r1} and {self.fastq_r2} are not of the same length!")

    def get_file_handles(self) -> List[TextIO]:
        """
        Opens the FASTQ files for reading and returns file handles

        :return: a list of file handles to the opened FASTQ file(s)
        """
        handles: List[TextIO] = [self._file_parser_handle(self.fastq_r1)]
        if not self.is_single_end:
            handles.append(self._file_parser_handle(self.fastq_r2))
        return handles

    @classmethod
    def open_func(cls, input_file) -> Callable:
        """
        Checks which function should be used to open an input_file, i.e. whether it is
        gzip-compressed or uncompressed

        :param input_file: the path to the file to be opened
        :return: a callable to be used for opening the input_file: this will be gzip.open for
        files having a .gz or .gzip file extension and the built-in open otherwise
        """
        gzip_extensions = [".gz", ".gzip"]
        return gzip.open if any([input_file.endswith(ext) for ext in gzip_extensions]) else open

    def _file_parser_handle(self, input_file: str) -> TextIO:
        with self.open_func(input_file)(input_file, 'rt') as file_handle:
            yield from file_handle


class SampleSheetParser:
    """
    Class that handles parsing of a simple tab-separated sample sheet
    """

    def __init__(self, samplesheet_file: str):
        """
        Create a new SampleSheetParser instance

        :param samplesheet_file: path to an existing tab-separated sample sheet
        """
        self.samplesheet_file: str = samplesheet_file

    def get_file_handle(self) -> TextIO:
        with open(self.samplesheet_file, 'rt') as fh:
            yield from fh

    def get_barcode_to_sample_mapping(self) -> Dict[str, str]:
        """
        Parse the sample sheet and return a barcode-to-sample mapping. For dual-index samples,
        the barcode will be created by concatenating the i5 and i7 barcodes using '+' as the
        delimiter

        :return: a dict having barcodes as keys and the corresponding sample-id as value
        """
        barcode_to_sample: Dict[str, str] = dict()
        for row in self.get_file_handle():
            pcs: List[str] = row.strip().split("\t")
            barcode_to_sample["+".join(pcs[1:])] = pcs[0]
        return barcode_to_sample
