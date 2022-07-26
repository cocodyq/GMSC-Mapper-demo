import subprocess
import argparse
import sys
import os
from os import path, makedirs
import pandas as pd
import tempfile

_ROOT = path.abspath(path.dirname(__file__))

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
                        required=False,
                        help='Path to the GMSC database file',
                        dest='database',
                        default=path.join(_ROOT, 'example/exampledb.dmnd'))

    parser.add_argument('--habitat', '--habitat',
                        required=False,
                        help='Path to the habitat file',
                        dest='habitat',
                        default=path.join(_ROOT, 'example/ref_habitat.txt'))
    parser.add_argument('--nohabitat','--nohabitat',action='store_true', help='Use this if no need to annotate habitat')

    parser.add_argument('--taxonomy', '--taxonomy',
                        required=False,
                        help='Path to the taxonomy file',
                        dest='taxonomy',
                        default=path.join(_ROOT, 'example/ref_taxonomy.txt'))
    parser.add_argument('--notaxonomy', '--notaxonomy',action='store_true', help='Use this if no need to annotate taxonomy')

    parser.add_argument('--quality', '--quality',
                        required=False,
                        help='Path to the quality file',
                        dest='quality',
                        default=path.join(_ROOT, 'example/ref_quality.txt'))
    parser.add_argument('--noquality', '--noquality',action='store_true', help='Use this if no need to annotate quality')

    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output directory (will be created if non-existent)',
                        dest='output',
                        default = None)	
					
    return parser.parse_args()

def validate_args(args):
    def expect_file(f):
        if not os.path.exists(f):
            sys.stderr.write(f"GMSC-mapper Error: Expected file '{f}' does not exist\n")
            sys.exit(1)
    
    if args.genome_fasta is None and args.aa_input is None and args.nt_input is None:
        sys.stderr.write("GMSC-mapper Error: At least one of --input or --aa-genes or --nt_genes is necessary\n")
        sys.stderr.exit(1)
    elif args.genome_fasta is not None and args.aa_input is None and args.nt_input is None:
        expect_file(args.genome_fasta)
    elif args.aa_input is not None and args.genome_fasta is None and args.nt_input is None:
        expect_file(args.aa_input)
    elif args.nt_input is not None and args.genome_fasta is None and args.aa_input is None:
        expect_file(args.nt_input)
    else:
        sys.stderr.write("GMSC-mapper Error: --input or --aa-genes or --nt_genes shouldn't be assigned at the same time\n")
        sys.stderr.exit(1)

    if args.database is None:
        sys.stderr.write("GMSC-mapper Error: GMSC databse is necessary\n")
        sys.stderr.exit(1)
    else:
        expect_file(args.database) 
    
    if args.habitat is not None:
        expect_file(args.habitat)

    if args.taxonomy is not None:
        expect_file(args.taxonomy)

    if args.quality is not None:
        expect_file(args.quality)

def predict_smorf(args):
    print('Start smORF prediction...')

    outdir = path.join(args.output,"predicted_smorf")

    subprocess.check_call([
        'macrel','get-smorfs',
        '--fasta',args.genome_fasta,
        '--output',outdir,
        '--cluster',
        '--keep-fasta-headers'])
    print('\nsmORF prediction has done.\n')
    return path.join(args.output,"predicted_smorf/macrel.out.smorfs.faa")

#should change parameter
def mapdb_diamond(queryfile,database,resultdir):
    print('Start smORF mapping...')

    #queryfile = path.join(args.output,"predicted_smorf/macrel.out.smorfs.faa")
    resultfile = path.join(resultdir,"diamond.out.smorfs.tsv")

    subprocess.check_call([
        'diamond','blastp',
        '-q',queryfile,
        '-d',database,
        '-o',resultfile,
        '--more-sensitive',
        '-e','0.00001',
        '--query-cover','90',
        '--subject-cover','90',
        '--outfmt','6','qseqid','full_qseq','qlen','sseqid','full_sseq','slen','pident','length','evalue','qcovhsp','scovhsp',
        '-p','64'])  
    print('\nsmORF mapping has done.\n')

def generate_fasta(queryfile,resultdir):
    print('Start smORF fasta file generating...')

    import pandas as pd
    from fasta import fasta_iter

    #queryfile = path.join(args.output,"predicted_smorf/macrel.out.smorfs.faa")
    resultfile = path.join(resultdir,"diamond.out.smorfs.tsv")
    fastafile = path.join(resultdir,"mapped.smorfs.faa")

    result = pd.read_csv(resultfile, sep='\t',header=None)
    smorf_id = set(result.iloc[:, 0].tolist())
    
    with open(fastafile,"wt") as f:
        for ID,seq in fasta_iter(queryfile):
            if ID in smorf_id:
                f.write(f'>{ID}\n{seq}\n')
    print('\nsmORF fasta file generating has done.\n')

def habitat(args):
    from map_habitat import smorf_habitat
    print('Start habitat annotation...')
    smorf_habitat(args)
    print('\nhabitat annotation has done.\n')

def taxonomy(args,tmpdirname):
    from map_taxonomy import deep_lca
    print('Start taxonomy annotation...')
    deep_lca(args,tmpdirname)
    print('\ntaxonomy annotation has done.\n')

def quality(args):
    from map_quality import smorf_quality
    print('Start quality annotation...')
    smorf_quality(args)
    print('\nquality annotation has done.\n')

def main(args=None):
    if args is None:
        args = sys.argv
    args = parse_args(args)
    validate_args(args)

    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    with tempfile.TemporaryDirectory() as tmpdirname:
        try:
            if args.genome_fasta:
                queryfile = predict_smorf(args)
            if args.aa_input:
                queryfile = args.aa_input
            mapdb_diamond(queryfile,args.database,args.output)
            generate_fasta(queryfile,args.output)
            if not args.nohabitat:
                habitat(args)
            if not args.notaxonomy:
                taxonomy(args,tmpdirname)
            if not args.noquality:
                quality(args)				
        except Exception as e:
            sys.stderr.write('GMGC-mapper Error: ')
            sys.stderr.write(str(e))
            sys.stderr.write('\n')
            sys.exit(1)		

if __name__ == '__main__':    
    main(sys.argv)
