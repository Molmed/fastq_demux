
import json
import os
import click

from fastq_demux.demux import FastqDemultiplexer, DemultiplexResults
from fastq_demux.parser import FastqFileParser, SampleSheetParser
from fastq_demux.writer import FastqFileWriterHandler


@click.command()
@click.option(
    "--R1",
    required=True,
    type=click.Path(
        exists=True,
        readable=True,
        resolve_path=True),
    help="FASTQ file with R1 sequences")
@click.option(
    "--R2",
    required=False,
    type=click.Path(
        exists=True,
        readable=True,
        resolve_path=True),
    help="FASTQ file with R2 sequences")
@click.option(
    "--samplesheet",
    required=True,
    type=click.Path(
        exists=True,
        readable=True,
        resolve_path=True),
    help="Sample sheet containing the barcode to sample-id mapping")
@click.option(
    "--prefix",
    required=False,
    type=click.STRING,
    default="",
    help="prefix to use for naming demultiplexed FASTQ files")
@click.option(
    "--unknown-barcode",
    required=False,
    type=click.STRING,
    default="Unknown",
    show_default=True,
    help="barcode to use for non-matching barcodes")
@click.option(
    "--outdir",
    required=False,
    type=click.Path(
        exists=True,
        resolve_path=True,
        dir_okay=True,
        writable=True),
    default=os.getcwd(),
    show_default="Current Working Directory",
    help="output directory where demultiplexed FASTQ files should be written")
@click.option(
    "--no-gzip-compression",
    flag_value=True,
    show_default=True,
    help="skip gzip-compression of output FASTQ files"
)
def demultiplex(r1, r2, samplesheet, prefix, unknown_barcode, outdir, no_gzip_compression) -> None:
    samplesheet_parser = SampleSheetParser(samplesheet_file=samplesheet)
    barcode_to_sample_mapping = samplesheet_parser.get_barcode_to_sample_mapping()
    barcode_to_sample_mapping[unknown_barcode] = "barcode"
    file_writer_handler = FastqFileWriterHandler(
        prefix=prefix,
        outdir=outdir,
        is_single_end=(r2 is None),
        no_gzip_compression=no_gzip_compression)
    file_writer_handler.fastq_file_writers_from_mapping(barcode_to_sample_mapping)

    fastq_parser = FastqFileParser(fastq_r1=r1, fastq_r2=r2)
    demultiplexer: FastqDemultiplexer = FastqDemultiplexer(
        fastq_parser=fastq_parser,
        fastq_writer=file_writer_handler,
        unknown_barcode=unknown_barcode)
    results: DemultiplexResults = demultiplexer.demultiplex()
    stats_file: str = os.path.join(
        outdir,
        f"{prefix}demux_Stats.json")
    with open(stats_file, 'wt') as fh:
        json.dump(
            results.stats_json(
                barcode_to_sample_mapping=samplesheet_parser.get_barcode_to_sample_mapping()),
            fh,
            indent=2)

    def _print_counts(header, counts):
        print("\t".join(header))
        print(
            "\n".join([
                "\t".join([
                    str(v) for v in count])
                for count in counts]))

    _print_counts(
        ["known_barcode", "count", "percent"],
        results.summarize_counts(results.known_barcodes))
    _print_counts(
        ["unknown_barcode", "count", "percent"],
        results.summarize_counts(results.unknown_barcodes, n_values=10))

