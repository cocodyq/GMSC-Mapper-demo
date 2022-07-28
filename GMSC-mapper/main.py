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
                        help='Path to the input genome contig sequence FASTA file.',
                        dest='genome_fasta',
                        default = None)

    #how to process nucleotide sequences
    parser.add_argument('--nt-genes', '--nt_genes',
                        required=False,
                        help='Path to the input nucleotide sequence FASTA file.',
                        dest='nt_input',
                        default=None)

    parser.add_argument('--aa-genes', '--aa_genes',
                        required=False,
                        help='Path to the input amino acid sequence FASTA file.',
                        dest='aa_input',
                        default=None)

    parser.add_argument('--tool', '--tool',
                        required=False,
						choices=['diamond', 'mmseqs'],
                        help='Sequence alignment tool(Diamond/MMseqs).',
                        dest='tool',
                        default='diamond')

    parser.add_argument('--db', '--db',
                        required=False,
                        help='Path to the GMSC database file.',
                        dest='database',
                        default=None)

    parser.add_argument('-s', '--sensitivity',
                        required=False,
                        help='Sensitivity.',
                        dest='sensitivity',
                        default=None)

    parser.add_argument('--id', '--id',
                        required=False,
                        help='Minimum identity to report an alignment(range 0.0-1.0).',
                        dest='identity',
                        default=0.0)

    parser.add_argument('--cov', '--cov',
                        required=False,
                        help='Minimum coverage to report an alignment(range 0.0-1.0).',
                        dest='coverage',
                        default=0.9)

    parser.add_argument('-e', '--evalue',
                        required=False,
                        help='Maximum e-value to report alignments(default=0.00001).',
                        dest='evalue',
                        default=0.00001)

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

    parser.add_argument('-t', '--threads',
                        required=False,
                        help='Number of CPU threads(default=3).',
                        dest='threads',
                        default=3)		

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

    if args.database is not None:
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
def mapdb_diamond(args,queryfile):
    print('Start smORF mapping...')

    resultfile = path.join(args.output,"diamond.out.smorfs.tsv")

    subprocess.check_call([
        'diamond','blastp',
        '-q',queryfile,
        '-d',args.database,
        '-o',resultfile,
        args.sensitivity,
        '-e',str(args.evalue),
        '--id',str(args.identity*100),
        '--query-cover',str(args.coverage*100),
        '--subject-cover',str(args.coverage*100),
        #'--outfmt','6 qseqid full_qseq qlen sseqid full_sseq slen pident length evalue qcovhsp scovhsp',can't work 
        # Value 6 may be followed by a space-separated list of these keywords: how to accept this parameter from commond
        '--outfmt','6','qseqid','full_qseq','qlen','sseqid','full_sseq','slen','pident','length','evalue','qcovhsp','scovhsp',
        '-p',str(args.threads)])  

    print('\nsmORF mapping has done.\n')
    return resultfile

def mapdb_mmseqs(args,queryfile,tmpdir):
    print('Start smORF mapping...')
    
    querydb = path.join(tmpdir,"query.db")
    resultdb = path.join(tmpdir,"result.db")
    tmp = path.join(tmpdir,"tmp","")
    resultfile = path.join(args.output,"mmseqs.out.smorfs.tsv")

    subprocess.check_call([
        'mmseqs','createdb',queryfile,querydb]) 

    subprocess.check_call([
        'mmseqs','search',
        querydb,
        args.database,
        resultdb,
		tmp,
		'-s',str(args.sensitivity),
        '-e',str(args.evalue),
        '--min-seq-id',str(args.identity),
        '-c',str(args.coverage),
        '--threads',str(args.threads)])  

    subprocess.check_call([
        'mmseqs','convertalis',
        querydb,
        args.database,
        resultdb,
        resultfile,
        '--format-output',"query,qseq,qlen,target,tseq,tlen,fident,alnlen,evalue,qcov,tcov"])		

    print('\nsmORF mapping has done.\n')
    return resultfile

def generate_fasta(args,queryfile,resultfile):
    print('Start smORF fasta file generating...')

    import pandas as pd
    from fasta import fasta_iter

    fastafile = path.join(args.output,"mapped.smorfs.faa")

    result = pd.read_csv(resultfile, sep='\t',header=None)
    smorf_id = set(result.iloc[:, 0].tolist())
    
    with open(fastafile,"wt") as f:
        for ID,seq in fasta_iter(queryfile):
            if ID in smorf_id:
                f.write(f'>{ID}\n{seq}\n')
    print('\nsmORF fasta file generating has done.\n')

def habitat(args,resultfile):
    from map_habitat import smorf_habitat
    print('Start habitat annotation...')
    smorf_habitat(args,resultfile)
    print('\nhabitat annotation has done.\n')

def taxonomy(args,resultfile,tmpdirname):
    from map_taxonomy import deep_lca
    print('Start taxonomy annotation...')
    deep_lca(args,resultfile,tmpdirname)
    print('\ntaxonomy annotation has done.\n')

def quality(args,resultfile):
    from map_quality import smorf_quality
    print('Start quality annotation...')
    smorf_quality(args,resultfile)
    print('\nquality annotation has done.\n')

def main(args=None):
    if args is None:
        args = sys.argv
    args = parse_args(args)

    if args.tool == 'diamond':
        if args.database is None:
            args.database = path.join(_ROOT, 'example/example_diamond_db.dmnd')	
        if args.sensitivity is None:
            args.sensitivity = '--more-sensitive'
        if args.sensitivity == '1':
            args.sensitivity = '--fast'
        if args.sensitivity == '2':
            args.sensitivity = '--mid-sensitive'
        if args.sensitivity == '3':
            args.sensitivity = '--default'
        if args.sensitivity == '4':
            args.sensitivity = '--sensitive'
        if args.sensitivity == '5':
            args.sensitivity = '--more-sensitive'
        if args.sensitivity == '6':
            args.sensitivity = 'very-sensitive'
        if args.sensitivity == '7':
            args.sensitivity = 'ultra-sensitive'
    if args.tool == 'mmseqs':
        if args.database is None:
            args.database = path.join(_ROOT, 'example/example_mmseqs_db')
        if args.sensitivity is None:
            args.sensitivity = 5.7

    validate_args(args)

    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    with tempfile.TemporaryDirectory() as tmpdirname:
        try:
            if args.genome_fasta:
                queryfile = predict_smorf(args)
            if args.aa_input:
                queryfile = args.aa_input

            if args.tool == 'diamond':
                resultfile = mapdb_diamond(args,queryfile)
            if args.tool == 'mmseqs':
                resultfile = mapdb_mmseqs(args,queryfile,tmpdirname)

            generate_fasta(args,queryfile,resultfile)

            if not args.nohabitat:
                habitat(args,resultfile)
            if not args.notaxonomy:
                taxonomy(args,resultfile,tmpdirname)
            if not args.noquality:
                quality(args,resultfile)				
        except Exception as e:
            sys.stderr.write('GMGC-mapper Error: ')
            sys.stderr.write(str(e))
            sys.stderr.write('\n')
            sys.exit(1)		

if __name__ == '__main__':    
    main(sys.argv)
