
import os

from typing import Dict, List, Optional, TextIO
from fastq_demux.parser import FastqFileParser


class FastqFileWriter:
    """
    Class handling writing of FASTQ reads to a file
    """

    def __init__(self, fastq_file: str):
        """
        Create a new FastqFileWriter instance and open the supplied file for writing. If the file
        already exists, it will be overwritten

        :param fastq_file: path to the FASTQ file that should be opended for writing
        """
        self.fastq_file: str = fastq_file
        self.fastq_file_handle: Optional[TextIO] = None
        self.open_handle()

    def open_handle(self) -> None:
        """
        Open the FASTQ file for writing, using the appropriate open function (i.e. gzip.open or
        open) and assign the file handle to the instance attribute "fastq_file_handle"
        """
        if self.fastq_file_handle is None:
            self.fastq_file_handle = FastqFileParser.open_func(self.fastq_file)(
                self.fastq_file, "wt")

    def close_handle(self) -> None:
        if self.fastq_file_handle is not None:
            self.fastq_file_handle.close()

    def write_record(self, record: List[str]) -> None:
        """
        Write the supplied FASTQ record to this instance's opened file handle

        :param record: a list of 4 strings representing a FASTQ read
        """
        self.fastq_file_handle.write("\n".join(record + [""]))


class FastqFileWriterHandler:
    """
    Class for handling a set of FastqFileWriters, where separate writers are used for each barcode
    and read.
    """

    def __init__(
            self,
            prefix: str,
            outdir: str,
            is_single_end: bool,
            no_gzip_compression: bool):
        """
        Create a new FastqFileWriterHandler instance.

        :param prefix: a string prefix to use in FASTQ file names
        :param outdir: the path to an existing output directory where FASTQ files should be written
        :param is_single_end: a boolean indicating whether single-end data should be written
        :param no_gzip_compression: a boolean indicating whether the output FASTQ files should be
        gzip-compressed
        """
        self.prefix: str = prefix
        self.outdir: str = outdir
        self.is_single_end: bool = is_single_end
        self.no_gzip_compression: bool = no_gzip_compression

        self.fastq_file_writers: Optional[Dict[str, List[FastqFileWriter]]] = None

    def fastq_file_name_from_sample_id(
            self, sample_id: str, barcode: str) -> List[str]:
        """
        Construct the name of an output FASTQ file based on the sample_id and barcode

        :param sample_id: the name of the sample whose reads will be written
        :param barcode: the barcode of the sample whose reads will be written
        :return: a list of paths to the FASTQ file(s) where the supplied sample's reads should be
        written
        """
        ext: str = "fastq.gz" if not self.no_gzip_compression else "fastq"
        return [
            os.path.join(
                self.outdir,
                f"{self.prefix}{sample_id}_{barcode.replace('+', '-')}_R{read_no}.{ext}"
            ) for read_no in range(1, 3 - int(self.is_single_end))]

    def fastq_file_writers_from_sample_id(
            self, sample_id: str, barcode: str) -> List[FastqFileWriter]:
        """
        Create FastqFileWriter instances to be used for writing reads for the supplied sample

        :param sample_id: the name of the sample whose reads will be written
        :param barcode: the barcode of the sample whose reads will be written
        :return: a list of FastqFileWriter instances to be used for writing reads for the supplied
        sample
        """
        return [
            FastqFileWriter(fastq_file)
            for fastq_file in self.fastq_file_name_from_sample_id(
                sample_id,
                barcode)]

    def fastq_file_writers_from_mapping(
            self,
            barcode_to_sample_mapping: Dict[str, str]) -> Dict[str, List[FastqFileWriter]]:
        """
        Create FastqFileWriter instances to be used for writing from a barcode-to-sample mapping.

        :param barcode_to_sample_mapping: a dict having barcodes as keys and the corresponding
        sample-id as value
        :return: a dict having barcodes as keys and a list of FastqFileWriter instances as value
        """
        if not self.fastq_file_writers:
            self.fastq_file_writers: Dict[str, List[FastqFileWriter]] = dict()
            for barcode, sample_id in barcode_to_sample_mapping.items():
                self.fastq_file_writers[barcode] = \
                    self.fastq_file_writers_from_sample_id(sample_id, barcode)

        return self.fastq_file_writers

    def write_fastq_record(self, barcode: str, fastq_record: List[List[str]]) -> None:
        """
        Write the supplied FASTQ read (or read pair) using the FastqFileWriter corresponding to the
        supplied barcode

        :param barcode: the barcode of the sample whose reads will be written
        :param fastq_record: a 1- or 2-element list of FASTQ reads for single-end or paired-end,
        respectively. Each read is a list of 4 strings
        :raises KeyError: if no FastqFileWriters exist for the supplied barcode
        """
        for read_no, record in enumerate(fastq_record):
            self.fastq_file_writers[barcode][read_no].write_record(record)
