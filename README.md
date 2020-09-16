# FASTQ demux

[![Build Status](https://travis-ci.org/Molmed/fastq_demux.svg?branch=master)](https://travis-ci.org/Molmed/fastq_demux)
[![codecov](https://codecov.io/gh/Molmed/fastq_demux/branch/master/graph/badge.svg)](https://codecov.io/gh/Molmed/fastq_demux)
[![Maintainability](https://api.codeclimate.com/v1/badges/e3809abe804678af9d79/maintainability)](https://codeclimate.com/github/Molmed/fastq_demux/maintainability)

A simple program to demultiplex Illumina FASTQ files based on barcodes in the FASTQ headers


- [Installation](#installation)
- [Basic usage](#basic-usage)
    - [Sample sheet](#sample-sheet)
- [Docker](#docker)


## Installation

```
git clone https://github.com/Molmed/fastq_demux.git
pip install -r requirements.txt
```
Alternatively, use a docker image (see [below](#docker)).


## Basic usage

For detailed usage instructions, run the script with the `--help` switch:
```
fastq_demux --help
```

To demultiplex a FASTQ file or a pair of FASTQ files based on the barcodes present in the FASTQ headers, supply a
file with forward reads (with `--R1`), reverse reads (with `--R2`, if paired-end) and a tab-separated sample sheet 
providing a barcode-to-sample mapping (with `--samplesheet`). 

```
fastq_demux --R1 tests/dual-index-short_Undetermined_S0_L001_R1_001.fastq.gz --samplesheet tests/samplesheet.tsv
```

### Sample sheet
The sample sheet should have two or three columns for single or dual index reads, respectively. The columns are 
`SampleID`, `P5-index` and `P7-index`.

Here is an example for a sample sheet accompanying a dual-index FASTQ file:
```
Sample_1    GGGGGGGG    AGATCTCG
Sample_2    GAAGATTT    TTTACTCT
Sample_3    GAAGATTT    AAAACGCC
```


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
