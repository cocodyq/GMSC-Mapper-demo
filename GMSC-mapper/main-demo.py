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
    parser.add_argument('--db', '--db',
                        required=True,
                        help='Path to the GMSC database file',
                        dest='database',
                        default=None)
    parser.add_argument('--habitat', '--habitat',
                        required=False,
                        help='Path to the habitat file',
                        dest='habitat',
                        default=None)
    parser.add_argument('--taxonomy', '--taxonomy',
                        required=False,
                        help='Path to the taxonomy file',
                        dest='taxonomy',
                        default=None)
    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output directory (will be created if non-existent)',
                        dest='output',
                        default = None)	
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
    expect_file(args.database) 

    if args.genome_fasta is None and args.aa_input is None and args.nt_input is None:
        sys.stderr.write("GMSC-mapper Error: At least one of --input or --aa-genes or --nt_genes is necessary\n")
        sys.stderr.exit(1)
    if args.database is None:
        sys.stderr.write("GMSC-mapper Error: GMSC databse is necessary\n")
        sys.stderr.exit(1)

def predict_smorf(args):
    outdir = path.join(args.output,"predicted_smorf")

    subprocess.check_call([
        'macrel','get-smorfs',
        '--fasta',args.genome_fasta,
        '--output',outdir,
        '--cluster',
        '--keep-fasta-headers'])

#should change parameter
def mapdb_diamond(args):
    queryfile = path.join(args.output,"predicted_smorf/macrel.out.smorfs.faa")
    resultfile = path.join(args.output,"diamond.out.smorfs.tsv")

    subprocess.check_call([
        'diamond','blastp',
        '-q',queryfile,
        '-d',args.database,
        '-o',resultfile,
        '--more-sensitive',
        '-e','0.00001',
        '--query-cover','90',
        '--subject-cover','90',
        '--outfmt','6','qseqid','full_qseq','qlen','sseqid','full_sseq','slen','pident','length','evalue','qcovhsp','scovhsp',
        '-p','64'])  

#should be optional
def generate_fasta(args):
    import pandas as pd
    from fasta import fasta_iter

    queryfile = path.join(args.output,"predicted_smorf/macrel.out.smorfs.faa")
    resultfile = path.join(args.output,"diamond.out.smorfs.tsv")
    fastafile = path.join(args.output,"mapped.smorfs.faa")

    result = pd.read_csv(resultfile, sep='\t',header=None)
    smorf_id = set(result.iloc[:, 0].tolist())
    
    with open(fastafile,"wt") as f:
        for ID,seq in fasta_iter(queryfile):
            if ID in smorf_id:
                f.write(f'>{ID}\n{seq}\n')

def habitat(args):
    from map_habitat import smorf_habitat

    smorf_habitat(args)
    
def taxonomy(args):
    from map_taxonomy import deep_lca

    deep_lca(args)

def main(args=None):
    if args is None:
        args = sys.argv
    args = parse_args(args)

    predict_smorf(args)
    mapdb_diamond(args)
    generate_fasta(args)
    habitat(args)
    taxonomy(args)

if __name__ == '__main__':    
    main(sys.argv)
