import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import fastq_demux.parser
import fastq_demux.demux
import fastq_demux.writer
