# rRNAFinder
rRNAFinder is a Python tool designed to identify rRNA genes in the genome.

## Installation
This package relies on NCBI BLAST and BLAT. Please ensure these are installed prior to installing rRNAFinder.

To install rRNAFinder using pip, use the following command:
```
pip install rRNAFinder
```

## Usage
```
usage: rRNAFinder [-h] [-s {plants,animals,fungi,bacteria}] [-l LOG_FILE] [-o OUTPUT_DIR] fasta_file

Identify rRNA genes in the genome

positional arguments:
  fasta_file            a fasta file

optional arguments:
  -h, --help            show this help message and exit
  -s {plants,animals,fungi,bacteria}, --species {plants,animals,fungi,bacteria}
  -l LOG_FILE, --log_file LOG_FILE
                        path for log file (default: None)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        output prefix (default: rRNA_finder)
```

Currently only plants genomes are supported