
from __future__ import annotations
import gzip
from typing import Callable, Dict, Iterator, List, Optional, TextIO, Tuple


class FastqFileParser:
    """
    Class that handles parsing of a single-end FASTQ file or a paired-end FASTQ file pair
    """

    def __init__(
            self,
            fastq_r1: str,
            fastq_r2: Optional[str] = None,
            fastq_i1: Optional[str] = None,
            fastq_i2: Optional[str] = None):
        """
        Create a new FastqFileParser instance

        :param fastq_r1: the path to the FASTQ file with forward reads (R1)
        :param fastq_r2: for paired-end, the path to the FASTQ file with reverse reads (R2)
        :param fastq_i1: the path to an optional FASTQ file with sequenced index i7 reads (I1).
                    If specified, will override any i7 index sequence in the FASTQ headers
        :param fastq_i2: the path to an optional FASTQ file with sequenced index i5 reads (I1).
                    If specified, will override any i5 index sequence in the FASTQ headers
        """
        self.fastq_r1: str = fastq_r1
        self.fastq_r2: Optional[str] = fastq_r2
        self.fastq_i1: Optional[str] = fastq_i1
        self.fastq_i2: Optional[str] = fastq_i2
        self.is_single_end: bool = self.fastq_r2 is None
        self.i7_in_header: bool = self.fastq_i1 is None
        self.i5_in_header: bool = self.fastq_i2 is None

    @classmethod
    def create_fastq_file_parser(
            cls,
            fastq_r1: str,
            fastq_r2: Optional[str] = None,
            fastq_i1: Optional[str] = None,
            fastq_i2: Optional[str] = None) -> FastqFileParser:
        if fastq_i1 is None:
            return FastqFileParserHeaderIndex(
                fastq_r1=fastq_r1,
                fastq_r2=fastq_r2
            )
        if fastq_i2 is None:
            return FastqFileParserI1(
                fastq_r1=fastq_r1,
                fastq_r2=fastq_r2,
                fastq_i1=fastq_i1
            )
        return FastqFileParserI2(
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            fastq_i1=fastq_i1,
            fastq_i2=fastq_i2
        )

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

    def fastq_records(self) -> Iterator[Tuple[List[List[str]], str]]:
        """
        Iterates over the FASTQ files

        :return: a generator which will yield FASTQ records. Each yielded element will be a list of
        length 1 for single-end and 2 for paired-end. Each list element will in turn be a list of
        4 strings, representing the 4 lines from the FASTQ file making up the FASTQ read
        :raises Exception: if, for paired-end, the supplied FASTQ files don't contain the same
        number of reads
        """
        read_handles, index_handles = self.get_file_handles()
        try:
            while True:
                yield self.records_from_handles(
                    read_handles=read_handles, index_handles=index_handles)
        except StopIteration:
            pass
        if not self.ensure_filehandles_empty(read_handles + index_handles):
            msg = f"All FASTQ files are not of the same length:\n"
            msg += "\n".join(
                [
                    f"\t{fq}"
                    for fq in (self.fastq_r1, self.fastq_r2, self.fastq_i1, self.fastq_i2)
                    if fq is not None])
            raise Exception(msg)

    def get_file_handles(self) -> Tuple[List[TextIO], List[TextIO]]:
        """
        Opens the FASTQ files for reading and returns file handles

        :return: a list of file handles to the opened FASTQ file(s)
        """
        handles: Tuple[List[TextIO], List[TextIO]] = ([self._file_parser_handle(self.fastq_r1)], [])
        for tuple_index, flag, fqfile in zip(
                [0, 1, 1],
                [self.is_single_end, self.i7_in_header, self.i5_in_header],
                [self.fastq_r2, self.fastq_i1, self.fastq_i2]):
            if not flag:
                handles[tuple_index].append(self._file_parser_handle(fqfile))

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

    @classmethod
    def barcode_from_record(cls, record: List[List[str]]) -> str:
        raise NotImplementedError

    def records_from_handles(
            self,
            read_handles: List[TextIO],
            index_handles: List[TextIO]) -> Tuple[List[List[str]], str]:
        return [self.record_from_handle(handle) for handle in read_handles], \
                self.barcode_from_record(
                    [self.record_from_handle(handle) for handle in index_handles])

    @classmethod
    def record_from_handle(cls, handle: TextIO) -> List[str]:
        return [next(handle).strip() for _ in range(4)]


class FastqFileParserHeaderIndex(FastqFileParser):

    @classmethod
    def barcode_from_record(cls, record: List[List[str]]) -> str:
        """
        Extract the barcode from a FASTQ read header. The header is expected to be ":"-delimited
        according to the bcl2fastq output format and the barcode is expected to be the last element
        when splitting on ":"

        :param record: a list of strings representing a set of FASTQ reads. Only the first element
        (the header string) of the first read will be accessed
        :return: the barcode string parsed from the header
        """
        return record[0][0].split(":")[-1]

    def records_from_handles(
            self,
            read_handles: List[TextIO],
            index_handles: List[TextIO]) -> Tuple[List[List[str]], str]:
        records = [self.record_from_handle(handle) for handle in read_handles]
        return records, self.barcode_from_record(records)


class FastqFileParserI1(FastqFileParser):

    @classmethod
    def barcode_from_record(cls, record: List[List[str]]) -> str:
        """
        Extract the barcode from a FASTQ index file. The single-index sequence will be in the last
        record in the list passed as input.

        :param record: a list of strings representing a set of FASTQ reads.
        :return: the barcode string present as a sequence in the FASTQ index file
        """
        return record[-1][1]


class FastqFileParserI2(FastqFileParser):

    @classmethod
    def barcode_from_record(cls, record: List[List[str]]) -> str:
        """
        Extract the barcode from a pair of FASTQ index files. The dual-index sequence will be in
        the two last records in the list passed as input.

        :param record: a list of strings representing a set of FASTQ reads.
        :return: the barcode string present as a sequence in the FASTQ index files
        """
        return f"{record[-2][1]}+{record[-1][1]}"


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
            pcs: List[str] = [r.strip() for r in row.strip().split("\t")]
            if len(pcs) > 1 and all([len(p) > 0 for p in pcs]):
                barcode_to_sample["+".join(pcs[1:])] = pcs[0]
        return barcode_to_sample
