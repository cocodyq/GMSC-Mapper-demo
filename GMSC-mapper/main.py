import subprocess
import os
import pandas as pd
    
def predict_smorf(contig,outdir):
    subprocess.check_call([
        'macrel','get-smorfs',
        '--fasta',contig,
        '--output',outdir,
        '--cluster',
        '--keep-fasta-headers'])

def map_diamond(outdir,outfile,dbfile):
    subprocess.check_call([
        'diamond','blastp',
        '-q',f'{outdir}/macrel.out.smorfs.faa',
        '-d',dbfile,
        '-o',f'{outdir}/{outfile}',
        '--more-sensitive',
        '-e','0.001',
        '--outfmt','6','qseqid','full_qseq','qlen','sseqid','full_sseq','slen','pident','length','evalue','qcovhsp','scovhsp',
        '-p','64'])  

def filter_cutoff(infile,outfile):
    result = pd.read_csv(infile, sep='\t',header=None)
    result = result[(result.iloc[:,6] > 90) & (result.iloc[:,8] < 0.00001)]
    result.to_csv(outfile,sep='\t',index=False,header=False)
    return outfile

def generate_fasta(infile,outfile):
    result = pd.read_csv(infile, sep='\t',header=None)
    result = result.drop_duplicates([0],keep='first')
    with open(outfile,"wt") as f:
        for index, row in result.iterrows():
            f.write(f'>{row[0]}\n{row[1]}\n')

def store_habitat(habitat_file):
    import lzma
    habitat_dict = {}
    with lzma.open(habitat_file,"rt") as f:
        for line in f:
            seq,habitat = line.strip().split("\t")
            habitat_dict[seq] = habitat
    return habitat_dict

def map_habitat(habitat_file,infile,outfile):
    habitat_dict = store_habitat(habitat_file)
    out = open(outfile,"wt")
    with open(infile,"rt") as f:
        for line in f:
            line = line.strip()
            linelist = line.split("\t")
            out.write(f'{line}\t{habitat_dict[linelist[3]]}\n')
    out.close()
    return outfile
    
def dedup_habitat(infile,outfile1):
    out = open(outfile1,"wt")
    query_habitat = {}
    with open(infile,"rt") as f:
        for line in f:
            linelist = line.strip().split("\t")
            if linelist[0] not in query_habitat.keys():
                query_habitat[linelist[0]] = set()
                single_habitat_list = linelist[11].split(",")
                for single_habitat in single_habitat_list:
                    query_habitat[linelist[0]].add(single_habitat)
            else:
                single_habitat_list = linelist[11].split(",")
                for single_habitat in single_habitat_list:
                    query_habitat[linelist[0]].add(single_habitat)
    for key,value in query_habitat.items():
        value = list(value)
        habitat = ",".join(value)
        out.write(f'{key}\t{habitat}\n')

def add_habitat(habitat_file,infile,outfile1,outfile2):
    outfile1 = map_habitat(habitat_file,infile,outfile1)
    dedup_habitat(outfile1,outfile2)
    return outfile1

def store_taxonomy(taxonomy_file):
    import lzma
    taxonomy_dict = {}
    with lzma.open(taxonomy_file,"rt") as f:
        for line in f:
            linelist = line.strip().split("\t",1)
            if len(linelist) == 2:
                taxonomy_dict[linelist[0]] = linelist[1]
            else:
                taxonomy_dict[linelist[0]] = ""
    return taxonomy_dict

def map_taxonomy(taxonomy_file,infile,outfile):
    taxonomy_dict = store_taxonomy(taxonomy_file)
    out = open(outfile,"wt")
    with open(infile,"rt") as f:
        for line in f:
            line = line.strip()
            linelist = line.split("\t")
            taxa = taxonomy_dict[linelist[3]].replace("\t",";")
            out.write(f'{line}\t{taxa}\n') 
    return outfile
            
def deep_lca(infile,outfile):
    cluster = {}
    taxa = {}
    change = {}
    lastrank = ""

    out = open(outfile, "wt")
    with open(infile,"rt") as f:
        for line in f:
            line = line.strip()
            linelist = line.split("\t")   
            if cluster:
                if linelist[0] in cluster.keys():
                    cluster[linelist[0]].append(linelist[3])
                else:
                    for key,value in cluster.items():
                        change[key] = []
                        for rank in range(7):
                            flag = 1
                            for i in range(len(value)):
                                # If a sequence has taxonomy
                                if taxa[value[i]]:
                                    #If a sequence has enough taxonomy to compare,if blank,continue
                                    if len(taxa[value[i]]) >= rank+1:
                                        if lastrank == "":
                                            lastrank = taxa[value[i]][rank]
                                        else:
                                            #If the current taxonomy is same as last one
                                            if taxa[value[i]][rank] != lastrank:
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
                                    change[key].append(lastrank)
                                    lastrank = ""
                            #for sequence doesn't have taxonomy
                                else:
                                    lastrank = ""
                                    break
                            else:
                                break  
                        #after all rank circle,output results
                        if change[key]:       
                            taxonomy = ";".join(change[key])
                            out.write(key+"\t"+taxonomy+"\n")
                        #for sequence doesn't have taxonomy
                        else:
                            out.write(key+"\n")                      
                    cluster = {}
                    taxa = {}
                    change = {}
                    cluster[linelist[0]] = [linelist[3]]
            else:
                cluster[linelist[0]] = [linelist[3]]
                
            taxa[linelist[3]] = []
            if len(linelist) > 12:
                for tax in linelist[12].split(";"):
                    taxa[linelist[3]].append(tax)
                
    for key,value in cluster.items():
        change[key] = []
        for rank in range(7):
            flag = 1
            for i in range(len(value)):
                if taxa[value[i]]:
                    if len(taxa[value[i]]) >= rank+1:
                        if lastrank == "":
                            lastrank = taxa[value[i]][rank]
                        else:
                            if taxa[value[i]][rank] != lastrank:
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
                    change[key].append(lastrank)
                    lastrank = ""
                else:
                    lastrank = ""
                    break   
            else:
                break   
        if change[key]:       
            taxonomy = ";".join(change[key])
            out.write(key+"\t"+taxonomy+"\n")
        else:
            out.write(key+"\n")                                  
    out.close()    
    
def add_taxonomy(taxonomy_file,infile,outfile1,outfile2):
    outfile1 = map_taxonomy(taxonomy_file,infile,outfile1)
    deep_lca(outfile1,outfile2)
    return outfile1

def map_quality(qulity_file,infile,outfile1):
    quality =  set()
    out = open(outfile1,"wt")
    with open(qulity_file,"rt") as f1:
        for line in f1:
            line = line.strip()
            quality.add(line)
    with open(infile,"rt") as f2:
        for line in f2:
            line = line.strip()
            linelist = line.split("\t")
            if len(linelist) == 13:
                if linelist[3] in quality:
                    out.write(f'{line}\thigh quality\n')
                else:
                    out.write(f'{line}\tlow quality\n')
            else:
                if linelist[3] in quality:
                    out.write(f'{line}\t\thigh quality\n')
                else:
                    out.write(f'{line}\t\tlow quality\n')                              
    out.close()
    
def cal_quality(infile,outfile2):
    quality = {}
    out = open(outfile2,"wt")
    with open(infile,"rt") as f1:
        for line in f1:
            linelist = line.strip().split("\t")
            if linelist[0] in quality.keys():
                quality[linelist[0]].add(linelist[13])   
            else:
                quality[linelist[0]] = set()
                quality[linelist[0]].add(linelist[13])
    for key,value in quality.items():
        if "high quality" in value:
            out.write(f'{key}\thigh quality\n')
        else:
            out.write(f'{key}\tlow quality\n')
            
def add_quality(qulity_file,infile,outfile1,outfile2):
    map_quality(qulity_file,infile,outfile1)
    cal_quality(outfile1,outfile2)

def main():
    RESULT_FILE = os.path.join(OUTPUT_PATH, OUTPUT_FILE)
    
    predict_smorf(INPUT_FILE,OUTPUT_PATH)
    
    map_diamond(OUTPUT_PATH,OUTPUT_FILE,DB_FILE)
    
    RESULT_FILE = filter_cutoff(RESULT_FILE,RESULT_FILE.replace(".tsv",".filtered.tmp"))
    
    generate_fasta(RESULT_FILE,RESULT_FILE.replace(".tmp",".faa"))
    
    RESULT_FILE = add_habitat(HABITAT_FILE,RESULT_FILE,RESULT_FILE.replace(".tmp",".habitat.tmp"),f'{OUTPUT_PATH}/habitat.tsv')
    
    RESULT_FILE = add_taxonomy(TAXONOMY_FILE,RESULT_FILE,RESULT_FILE.replace(".tmp",".taxa.tmp"),f'{OUTPUT_PATH}/taxonomy.tsv')
    
    add_quality(QUALITY_FILE,RESULT_FILE,RESULT_FILE.replace(".tmp",".quality.tsv"),f'{OUTPUT_PATH}/quality.tsv')
    
    
INPUT_FILE = "~/GMSC/mapper/wilmes/Wilmes_Shared_smORFs/CAMPI/Fecal.contigs.fa"
OUTPUT_PATH = "~/GMSC/mapper/wilmes/final_results/old/sihumix"
OUTPUT_FILE = "diamond.result.tsv"
DB_FILE = "~/GMSC/mapper/diamond/index/90AA_GMSC"   
HABITAT_FILE = "~/GMSC/habitat/multi/id90/90AA_ref_multi_general_habitat.tsv.xz"
TAXONOMY_FILE = "~/GMSC/data/frozen/90AA_taxonomy.tsv.xz"
QUALITY_FILE = "~/GMSC/quality/allpass_90_delcoor.txt"

if __name__ == '__main__':    
    main()
