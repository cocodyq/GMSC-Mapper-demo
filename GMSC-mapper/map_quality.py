from os import path
import pandas as pd

def smorf_quality(args):
    result_file = path.join(args.output,"diamond.out.smorfs.tsv")
    quality_file = path.join(args.output,"quality.out.smorfs.tsv")	

    result = pd.read_csv(result_file, sep='\t',header=None)
    result.columns = ['qseqid','full_qseq','qlen','sseqid','full_sseq','slen','pident','length','evalue','qcovhsp','scovhsp']
    ref_quality = pd.read_csv(args.quality, sep='\t',header=None)
    ref_quality.columns = ['sseqid']
    ref_quality['quality'] = 'high quality'

    output = pd.merge(result,ref_quality,how='left')[['qseqid', 'quality']].fillna('low quality')
    output = output.groupby('qseqid',as_index=False,sort=False).agg({'quality':lambda x : ','.join(x.drop_duplicates())})
    rule = {"high quality":0,"low quality":1}
    output['quality'] = output['quality'].apply(lambda x: sorted(x.split(','),key=lambda x:rule[x])[0])
    output.to_csv(quality_file,sep='\t',index=False)
    return quality_file