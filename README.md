# FASTQ demux

[![Build Status](https://travis-ci.org/Molmed/fastq_demux.svg?branch=master)](https://travis-ci.org/Molmed/fastq_demux)
[![codecov](https://codecov.io/gh/Molmed/fastq_demux/branch/master/graph/badge.svg)](https://codecov.io/gh/Molmed/fastq_demux)
[![Maintainability](https://api.codeclimate.com/v1/badges/e3809abe804678af9d79/maintainability)](https://codeclimate.com/github/Molmed/fastq_demux/maintainability)

A simple program to demultiplex Illumina FASTQ files based on barcodes in the FASTQ headers

## Docker
To build a docker image:

```
git clone https://github.com/Molmed/fastq_demux.git
docker build -t fastq_demux:master .
```

Run the docker image without arguments to see usage:

```
docker run fastq_demux:master
```

Example usage:

```
docker run fastq_demux:master \
--R1 /code/tests/dual-index-short_Undetermined_S0_L001_R1_001.fastq.gz \
--samplesheet /code/tests/samplesheet.tsv
```

Another example using a data directory mounted into the container from the local filesystem:

```
docker run -v $(pwd)/tests:/data fastq_demux:master \
--R1 /data/dual-index-short_Undetermined_S0_L001_R1_001.fastq.gz \
--samplesheet /data/samplesheet.tsv
```
