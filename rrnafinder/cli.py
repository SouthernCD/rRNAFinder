import argparse
from rrnafinder.pipeline import rRNAFinder_main

def main():

    # argument parse
    parser = argparse.ArgumentParser(
        prog='rRNAFinder', description='Find rRNA unit and all rRNA in genome\n'
    )


    parser.add_argument('fasta_file', type=str, help='a fasta file')
    parser.add_argument('-s', '--species', type=str, choices=['plants', 'animals', 'fungi', 'bacteria'], default='plants')
    parser.add_argument('-l', '--log_file', type=str,
                          help='path for log file (default: None)', default=None)
    parser.add_argument('-o', '--output_dir', type=str, help='output prefix (default: rRNA_finder)',
                          default='rRNA_finder')

    args = parser.parse_args()

    rRNAFinder_main(args)

if __name__ == '__main__':
    main()