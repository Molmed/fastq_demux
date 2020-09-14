
import click
import os
from collections import Counter

import fastq_demux.demux
from fastq_demux.parser import FastqFileParser, SampleSheetParser


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
def demultiplex(r1, r2, samplesheet, prefix, unknown_barcode, outdir):
    demultiplexer = fastq_demux.demux.Demultiplexer(
        fastq_file_parser=FastqFileParser(fastq_r1=r1, fastq_r2=r2),
        samplesheet_parser=SampleSheetParser(samplesheet_file=samplesheet),
        prefix=prefix,
        unknown_barcode=unknown_barcode,
        outdir=outdir)
    (known_counts, unknown_counts) = demultiplexer.demultiplex()
    print("\t".join(["known_barcode", "count", "percent"]))
    print(
        "\n".join([
            "\t".join([
                str(v) for v in counts])
            for counts in demultiplexer.format_counts(known_counts)]))

    print("\t".join(["unknown_barcode", "count", "percent"]))
    print(
        "\n".join([
            "\t".join([
                str(v) for v in counts])
            for counts in demultiplexer.format_counts(unknown_counts)]))


if __name__ == '__main__':
    demultiplex()
