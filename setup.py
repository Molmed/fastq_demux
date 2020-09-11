from setuptools import setup, find_packages
from fastq_demux import __version__

setup(
    name='fastq_demux',
    version=__version__,
    description="A simple program to demultiplex Illumina FASTQ files based on barcodes in the fastq headers",
    long_description="A simple program to demultiplex Illumina FASTQ files based on barcodes in the fastq headers",
    keywords=['bioinformatics', 'illumina', 'demultiplexing', 'FASTQ'],
    author='Pontus Larsson, SNP&SEQ Technology Platform, Uppsala University',
    author_email='pontus.larsson@medsci.uu.se',
    url="https://www.github.com/",
    download_url='https://github.com/Molmed/checkQC/archive/{}.tar.gz'.format(__version__),
    install_requires=[
        "click"],
    packages=find_packages(exclude=["tests*"]),
    test_suite="tests",
    package_data={},
    include_package_data=True,
    license='GPLv3',
    entry_points={
        'console_scripts': ['fastq_demux = fastq_demux.fastq_demux:demultiplex']
    },
)
