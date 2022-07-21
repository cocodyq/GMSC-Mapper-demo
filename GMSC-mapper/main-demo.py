import subprocess
import argparse
import sys
import os
from os import path, makedirs
import pandas as pd
    
def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='GMSC-mapper')
    parser.add_argument('-i', '--input',
                        required=False,
                        help='Path to the input genome FASTA file.',
                        dest='genome_fasta',
                        default = None)
    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output directory (will be created if non-existent)',
                        dest='output',
                        default = None)
    parser.add_argument('--nt-genes', '--nt_genes',
                        required=False,
                        help='Path to the input DNA gene file (FASTA format)',
                        dest='nt_input',
                        default=None)
    parser.add_argument('--aa-genes', '--aa_genes',
                        required=False,
                        help='Path to the input amino acid gene file (FASTA format)',
                        dest='aa_input',
                        default=None)
    return parser.parse_args()

def validate_args(args):
    def expect_file(f):
        if f is not None:
            if not os.path.exists(f):
                sys.stderr.write(f"GMSC-mapper Error: Expected file '{f}' does not exist\n")
                sys.exit(1)
    expect_file(args.genome_fasta)
    expect_file(args.aa_input)
    expect_file(args.nt_input)
    if args.genome_fasta is None and args.aa_input is None and args.nt_input:
        sys.stderr.write("GMSC-mapper Error: At least one of --input or --aa-genes or --nt_genes is necessary\n")
        sys.stderr.exit(1)

def predict_smorf(args):
	outdir = path.join(args.output,"predicted_smorf")
	subprocess.check_call([
        'macrel','get-smorfs',
        '--fasta',args.genome_fasta,
        '--output',outdir,
        '--cluster',
        '--keep-fasta-headers'])

def main(args=None):
	if args is None:
		args = sys.argv
	args = parse_args(args)
	predict_smorf(args)

if __name__ == '__main__':    
    main(sys.argv)