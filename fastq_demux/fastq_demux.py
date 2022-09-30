
import json
import os
import click

from fastq_demux.demux import FastqDemultiplexer, DemultiplexResults
from fastq_demux.parser import FastqFileParser, SampleSheetParser
from fastq_demux.writer import FastqFileWriterHandler


def _validate_index_files(ctx, param, value):
    if param.name == 'i2' and 'i1' not in ctx.params:
        raise click.UsageError(
            "you must also specify an I1 index file specifying an I2 index file. Note that the "
            "order of the options matter: --I1 should occur before --I2")


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
    "--I1",
    required=False,
    type=click.Path(
        exists=True,
        readable=True,
        resolve_path=True),
    help="FASTQ file with I1 (i7) sequences. " \
         "If specified, will replace any index present in the FASTQ headers")
@click.option(
    "--I2",
    required=False,
    type=click.Path(
        exists=True,
        readable=True,
        resolve_path=True),
    callback=_validate_index_files,
    help="FASTQ file with I2 (i5) sequences. " \
         "If specified, will replace any index present in the FASTQ headers")
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
def demultiplex(r1, r2, i1, i2, samplesheet, prefix, unknown_barcode, outdir, no_gzip_compression) -> None:
    samplesheet_parser = SampleSheetParser(samplesheet_file=samplesheet)
    barcode_to_sample_mapping = samplesheet_parser.get_barcode_to_sample_mapping()
    barcode_to_sample_mapping[unknown_barcode] = "barcode"
    file_writer_handler = FastqFileWriterHandler(
        prefix=prefix,
        outdir=outdir,
        is_single_end=(r2 is None),
        no_gzip_compression=no_gzip_compression)
    file_writer_handler.fastq_file_writers_from_mapping(barcode_to_sample_mapping)

    fastq_parser = FastqFileParser.create_fastq_file_parser(
        fastq_r1=r1, fastq_r2=r2, fastq_i1=i1, fastq_i2=i2)
    results = DemultiplexResults(barcode_to_sample_mapping=barcode_to_sample_mapping)
    demultiplexer: FastqDemultiplexer = FastqDemultiplexer(
        fastq_parser=fastq_parser,
        fastq_writer=file_writer_handler,
        demultiplex_results=results,
        unknown_barcode=unknown_barcode)
    demultiplexer.demultiplex()
    stats_file: str = os.path.join(
        outdir,
        f"{prefix}demux_Stats.json")
    with open(stats_file, 'wt') as fh:
        json.dump(
            results.stats_json(),
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
        results.summarize_counts(barcode_set="known"))
    _print_counts(
        ["unknown_barcode", "count", "percent"],
        results.summarize_counts(barcode_set="unknown", n_values=10))
