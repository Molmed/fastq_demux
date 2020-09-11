
import click
import os
from collections import Counter

import fastq_demux.demux


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
        fastq_file_r1=r1,
        fastq_file_r2=r2,
        samplesheet=samplesheet,
        prefix=prefix,
        unknown_barcode=unknown_barcode,
        outdir=outdir
    )
    (known_counts, unknown_counts) = demultiplexer.demultiplex()
    print("\t".join(["known_barcode", "count", "percent"]))
    print("\n".join(["\t".join([barcode, str(count), f"{round(100.*count/sum(known_counts.values()), 1)}%"]) for barcode, count in known_counts.most_common()]))
    print("\t".join(["unknown_barcode", "count", "percent"]))
    print("\n".join(["\t".join([barcode, str(count), f"{round(100.*count/sum(unknown_counts.values()), 1)}%"]) for barcode, count in unknown_counts.most_common(10)]))


if __name__ == '__main__':
    demultiplex()
