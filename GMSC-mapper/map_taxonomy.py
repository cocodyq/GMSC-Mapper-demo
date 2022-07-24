# how to solve memory/time consuming
# adapt to different format(.xz .gz)
from os import path
import pandas as pd

def smorf_taxonomy(args):
    result_file = path.join(args.output,"diamond.out.smorfs.tsv")
    taxonomy_file = path.join(args.output,"taxonomy.out.smorfs.tmp.tsv")	

    result = pd.read_csv(result_file, sep='\t',header=None)
    result.columns = ['qseqid','full_qseq','qlen','sseqid','full_sseq','slen','pident','length','evalue','qcovhsp','scovhsp']
    ref_taxonomy = pd.read_csv(args.taxonomy, sep='\t',header=None)
    ref_taxonomy.columns = ['sseqid','taxonomy']

    output = pd.merge(result,ref_taxonomy,how='left')[['qseqid', 'sseqid','taxonomy']]
    output.to_csv(taxonomy_file,sep='\t',index=False)
    return taxonomy_file

def deep_lca(args):
    taxonomy_file = smorf_taxonomy(args)
    taxonomy_dlca_file = path.join(args.output,"taxonomy.out.smorfs.dlca.tsv")
	
    cluster = {}
    taxa = {}
    change = {}
    lastrank = ""

    out = open(taxonomy_dlca_file, "wt")
    with open(taxonomy_file,"rt") as f:
        for line in f:
            linelist = line.strip().split("\t")
            qseqid = linelist[0]
            sseqid = linelist[1] 	
            if len(linelist) == 3:
                original_tax = linelist[2]

            if cluster:
                if qseqid in cluster.keys():
                    cluster[qseqid].append(sseqid)
                else:
                    for q_seqid,s_seqid_list in cluster.items():
                        change[q_seqid] = []
                        for rank in range(7):
                            flag = 1
                            for s_seqid in s_seqid_list:
                                # If a sequence has taxonomy
                                if taxa[s_seqid]:
                                    #If a sequence has enough taxonomy to compare,if blank,continue
                                    if len(taxa[s_seqid]) >= rank+1:
                                        if lastrank == "":
                                            lastrank = taxa[s_seqid][rank]
                                        else:
                                            #If the current taxonomy is same as last one
                                            if taxa[s_seqid][rank] != lastrank:
                                                flag = 0
                                                lastrank = ""
                                                break
                                            else:
                                                continue
                                    else:   
                                        continue
                                else:  
                                    continue
                            #add same taxonomy every rank
                            if flag == 1:
                                if lastrank != "":
                                    change[q_seqid].append(lastrank)
                                    lastrank = ""
                            #for sequence doesn't have taxonomy
                                else:
                                    lastrank = ""
                                    break
                            else:
                                break  
                        #after all rank circle,output results
                        if change[q_seqid]:       
                            taxonomy = ";".join(change[q_seqid])
                            out.write(q_seqid+"\t"+taxonomy+"\n")
                        #for sequence doesn't have taxonomy
                        else:
                            out.write(q_seqid+"\n")                      
                    cluster = {}
                    taxa = {}
                    change = {}
                    cluster[qseqid] = [sseqid]
            else:
                cluster[qseqid] = [sseqid]
                
            taxa[sseqid] = []
            if len(linelist) == 3:
                for tax in original_tax.split(";"):
                    taxa[sseqid].append(tax)
                
    for q_seqid,s_seqid_list in cluster.items():
        change[q_seqid] = []
        for rank in range(7):
            flag = 1
            for s_seqid in s_seqid_list:
                if taxa[s_seqid]:
                    if len(taxa[s_seqid]) >= rank+1:
                        if lastrank == "":
                            lastrank = taxa[s_seqid][rank]
                        else:
                            if taxa[s_seqid][rank] != lastrank:
                                flag = 0
                                lastrank = ""
                                break
                            else:
                                continue
                    else:
                        continue
                else:
                    continue
            if flag == 1:
                if lastrank != "":
                    change[q_seqid].append(lastrank)
                    lastrank = ""
                else:
                    lastrank = ""
                    break   
            else:
                break   
        if change[q_seqid]:       
            taxonomy = ";".join(change[q_seqid])
            out.write(q_seqid+"\t"+taxonomy+"\n")
        else:
            out.write(q_seqid+"\n")                                  
    out.close()  

'''
def store_taxonomy(args):
    import lzma
    taxonomy_dict = {}
    with open(args.taxonomy,"rt") as f:
        for line in f:
            linelist = line.strip().split("\t")
            if len(linelist) == 2:
                taxonomy_dict[linelist[0]] = linelist[1]
            else:
                taxonomy_dict[linelist[0]] = ""
    return taxonomy_dict

def map_taxonomy(args):
    taxonomy_dict = store_taxonomy(args)
    result_file = path.join(args.output,"diamond.out.smorfs.tsv")
    taxonomy_mapped_file = path.join(args.output,"taxonomy.out.smorfs.tsv")

    with open(taxonomy_mapped_file,"wt") as out:
        with open(result_file,"rt") as f:
            for line in f:
                line = line.strip()
                linelist = line.split("\t")
                taxa = taxonomy_dict[linelist[3]]
                out.write(f'{line}\t{taxa}\n') 
'''			