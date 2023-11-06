# rRNAFinder
A python tool for identify rRNA gene in plant genome.

## Installation
This package depended on ncbi blast, you should install it first.

Install rRNAFinder with pip:
```
pip install rRNAFinder
```

## Usage
```
usage: rRNAFinder [-h] [-s QUERY_RRNA_SEQ] [-b QUERY_RRNA_BED] [-l LOG_FILE] [-o OUTPUT_DIR] fasta_file

Find rRNA unit and all rRNA in genome

positional arguments:
  fasta_file            a fasta file

optional arguments:
  -h, --help            show this help message and exit
  -s QUERY_RRNA_SEQ, --query_rRNA_seq QUERY_RRNA_SEQ
                        a known rRNA unit seq from close species
  -b QUERY_RRNA_BED, --query_rRNA_bed QUERY_RRNA_BED
                        bed file for close rRNA unit seq
  -l LOG_FILE, --log_file LOG_FILE
                        path for log file (default: None)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        output prefix (default: rRNA_finder)
```