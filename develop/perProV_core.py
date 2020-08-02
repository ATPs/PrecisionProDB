import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import gzip
from io import StringIO
import numpy as np
from multiprocessing import Pool
from collections import Counter
import os
import time

#file_genome = '/projectsp/f_jx76_1/xiaolong/genome/human/GRCh38/GCA_000001405.28_GRCh38.p13_genomic.fna.gz'
#file_proteins = '/projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release33/gencode.v33.pc_translations.fa.gz'
#file_gtf = '/projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release33/gencode.v33.chr_patch_hapl_scaff.annotation.gtf.gz'
#file_mutations = '/projectsp/f_jx76_1/xiaolong/2020humanRefPr/GENCODE/gnomAD3AF0.01ExonEthnics/adj.csv.gz'
#outprefix = '/projectsp/f_jx76_1/xiaolong/2020humanRefPr/GENCODE/20200216EthnicProteins/all'
#cpu_counts = os.cpu_count()


def parse_mutations(file_mutations):
    '''
    return a dataframe of mutations
    '''
    # read in mutations
    df_mutations = pd.read_csv(file_mutations, sep='\t',low_memory=False)
    
    # change chr column to the format of single letter, not begin with 'chr'
    df_mutations['chr'] = df_mutations['chr'].astype(str)
    df_mutations['chr'] = df_mutations['chr'].apply(lambda x:x[3:] if x.startswith('chr') else x)
    return df_mutations

def parse_gtf_gencode(file_gtf):
    '''
    given a gtf file of gencode, filter only coding transcripts. 
    For the attribute field, convert them to columns of df_gtf
    '''
    # read in gtf file
    df_gtf = pd.read_csv(file_gtf, sep='\t',header=None, comment='#')
    columns_gtf = ['seqname','source','feature','start','end','score','strand','frame','attribute']
    df_gtf.columns = columns_gtf
    # change seqname to single letter if possible
    df_gtf['seqname'] = df_gtf['seqname'].astype(str)
    df_gtf['seqname'] = df_gtf['seqname'].apply(lambda x:x[3:] if x.startswith('chr') else x)
    ## expand annotation in "attribute" column
    tls = []
    for i in df_gtf['attribute']:
        tdc = {}
        es = i.strip(';').split(';')
        for e in es:
            k,v = e.split()
            v = v.strip('"')
            tdc[k] = v
        tls.append(tdc)
    tdf = pd.DataFrame(tls)
    df_gtf = pd.concat([df_gtf, tdf], axis=1)
    df_gtf = df_gtf.drop('attribute', axis=1)
    ## only keep those with protein sequences
    df_gtf = df_gtf[df_gtf['protein_id'].notna()].copy()

    return df_gtf

def parse_proteins_gencode(file_proteins):
    '''
    return a dataframe with file_proteins from gencode
    '''
    # get protein sequences
    if file_proteins.endswith('.gz'):
        tls =list(SeqIO.parse(gzip.open(file_proteins,'rt'), 'fasta'))
    else:
        tls =list(SeqIO.parse(open(file_proteins,'r'), 'fasta'))
    df_protein = pd.DataFrame([seq.id.split('|') + [str(seq.seq)] for seq in tls])
    df_protein.columns = ['protein_id', 'transcript_id', 'gene_id', 'havana_gene', 'havana_transcript', 'transcript_name', 'gene_name', 'AA_length', 'AA_seq']
    print('from the protein file, totally', len(tls), 'protein sequences. number of unique gene_id, transcript_id and protein_id are', df_protein.gene_id.unique().shape[0], df_protein.transcript_id.unique().shape[0], df_protein.protein_id.unique().shape[0])
    return df_protein

def parse_genome(file_genome):
    '''
    file_genome can be a gzip file
    given a genome file, return a dictionary. with sequence_id (1-22, XYM, then fragment names) as keys, a string of bases in upper case as values
    '''
    chromosomes = [str(i) for i in range(1,23)] + list('XYM')
    # read in genome to dictionary and keep only chromosomes
    if file_genome.endswith('.gz'):
        tdc = SeqIO.to_dict(SeqIO.parse(gzip.open(file_genome,'rt'), 'fasta'))
    else:
        tdc = SeqIO.to_dict(SeqIO.parse(open(file_genome,'r'), 'fasta'))
    dc_chrom = {}
    for k,v in tdc.items():
        for chromosome in chromosomes:
            if f'chromosome {chromosome}, GRCh38 reference primary assembly' in v.description:
                dc_chrom[chromosome] = str(v.seq).upper()
                break
        else:
            if f'Homo sapiens mitochondrion, complete genome' in v.description:
                dc_chrom['M'] = str(v.seq).upper()
            else:
                dc_chrom[k] = str(v.seq).upper()
    del tdc
    print('number of genomic fragments is :', len(dc_chrom))
    print('major chromosomes', [e for e in chromosomes if e in dc_chrom])
    return dc_chrom



def get_df_transcripts_from_df_gtf(df_gtf):
    # create dataframe of transcripts based on gtf file
    tdf = df_gtf[df_gtf['feature'] == 'CDS']
    df_transcript2 = tdf.drop_duplicates(subset='transcript_id', keep='first') # the same size as df_transcript
    df_transcript2 = df_transcript2[['seqname','gene_id', 'transcript_id', 'protein_id','gene_name','hgnc_id','strand']]
    print('from the gtf file,', df_transcript2.shape[0], 'transcripts can be translated.')
    return df_transcript2.copy()

def addAA_seq(df_transcript2, file_proteins):
    '''
    file_proteins is the file with annotated protein sequences
    add column AA_seq to df_transcript2, as the annotated protein sequences
    '''
    if file_proteins is None:
        df_transcript2['AA_seq'] = None
        return df_transcript2
    
    # get df_protein
    df_protein = parse_proteins_gencode(file_proteins)
    ## append protein sequence to df_transcript2
    tdf = df_protein[['transcript_id','AA_seq']]# cannot use protein_id, as some protein_id from same transcripts, not one to one due to alternative fragments in the genome
    print('number of unique transcripts in gtf file and protein sequence file are', df_transcript2.shape[0], tdf.shape[0])
    if tdf['transcript_id'].unique().shape[0] != tdf.shape[0]:
        print('something wrong. the transcript id in protein sequence file is not unique')
    df_transcript2 = df_transcript2.merge(tdf, left_on = 'transcript_id', right_on='transcript_id', how='left')#110809
    return df_transcript2


def getCDSplus(tdf):
    '''
    dc_chrom is the dict of genomic DNA
    given a transcript_id, get the DNA_seq based on the gtf annotation
    dcgtf_transcript is a dict with transcript_id as key, and gtf dataframe for that transcript as value
    '''
#    tdf = dcgtf_transcript[transcript_id]
    transcript_id, tdf = tdf
    tdf = tdf[tdf['feature'].isin(['exon','CDS'])]
    strand = tdf['strand'].iloc[0]
    if strand =='+':
        tdf = tdf.sort_values(by='start')# sort by location from small to large if positive strand
    else:
        tdf = tdf.sort_values(by='start', ascending=False) #negative strand, from large to small
    # remove rows before the first exon till the first CDS, then remove all CDS
    keep = []
    for i in range(tdf.shape[0]):
        if tdf.iloc[i]['feature'] != 'CDS':
            keep.append(False)
        else:
            keep.append(True)
            frame = tdf.iloc[i]['frame']
            break
    for j in range(i+1, tdf.shape[0]):
        if tdf.iloc[j]['feature'] == 'CDS':
            keep.append(False)
        else:
            if strand == '+':
                if tdf.iloc[j]['start'] >= tdf.iloc[i]['end']:
                    keep.append(True)
                else:
                    keep.append(False)
            else:
                if tdf.iloc[j]['end'] <= tdf.iloc[i]['start']:
                    keep.append(True)
                else:
                    keep.append(False)
    tdf = tdf[keep]
    tlocs = [(i-1,j) for i,j in zip(tdf['start'], tdf['end'])]
    tlocs = sorted(tlocs, key=lambda x:x[0]) # sort by location from small to large
    return transcript_id, frame, tlocs
    
def addCDSplus_notes(df_transcript2, df_gtf, file_genome, cpu_counts):
    '''
    get CDSplus sequences, frame, genomicLocs, genomicStart, genomicEnd for df_transcript2
    df_gtf is the gtf dataframe
    file_genome is genomic DNA sequence
    '''
    # get CDSplus sequences
    df_gtf_group = df_gtf.groupby('transcript_id')
    
    pool = Pool(cpu_counts)
    results = pool.map(getCDSplus, df_gtf_group)
    pool.close()
    pool.join()
    
    tdf = pd.DataFrame(results)
    tdf.columns = ['transcript_id','frame', 'genomicLocs']
    df_transcript2 = df_transcript2.merge(tdf, left_on = 'transcript_id', right_on='transcript_id', how='left')
    print('finishing get locs and frame')
    
    # get genome DNA sequences
    dc_chrom = parse_genome(file_genome)
    print('finish reading genome file')
    
    df_transcript2['CDSplus'] = df_transcript2.apply(lambda x:''.join([dc_chrom[x['seqname']][i:j] for i,j in x['genomicLocs']]), axis=1)
    df_transcript2['genomicStart'] = df_transcript2['genomicLocs'].apply(lambda x:x[0][0])
    df_transcript2['genomicEnd'] = df_transcript2['genomicLocs'].apply(lambda x:x[-1][1])
    return df_transcript2

def get_df_transcript2(file_gtf, file_proteins, file_genome, cpu_counts):
    '''
    get df_transcripts, add AA_seq, add CDSplus and other notes
    '''
    
    # get df_gtf
    df_gtf = parse_gtf_gencode(file_gtf)
    
    # get df_transcript2, with only coding genes
    df_transcript2 = get_df_transcripts_from_df_gtf(df_gtf)
    
    # add AA_seq, annotated protein sequence to df_transcript2
    df_transcript2 = addAA_seq(df_transcript2, file_proteins)
    
    # get CDSplus sequences, frame, genomicLocs, genomicStart, genomicEnd for df_transcript2
    df_transcript2 = addCDSplus_notes(df_transcript2, df_gtf, file_genome, cpu_counts)
    # set "transcript_id" as index
    df_transcript2 = df_transcript2.set_index('transcript_id')
    return df_transcript2

# assign mutations to each transcript. 
# add column 'mutations' to df_transcript2, which stores index on the mutations in df_mutations
def getMutations(transcript_id):
    '''
    given a row in df_transcript2 of transcript_id, return the index of mutations where it locates inside genomicStart and genomicEnd
    Note the genomicStart should be add 1, as it is index starting from 0
    '''
    r = df_transcript2.loc[transcript_id]
    genomicStart = r['genomicStart'] + 1
    genomicEnd = r['genomicEnd']
    chromosome = r['seqname']
    tdf = df_mutations[(df_mutations['chr'] == chromosome)]
    tdf = tdf[(tdf['pos'] >= genomicStart) & (tdf['pos'] <= genomicEnd)]
    return list(tdf.index)


def translateCDSplus(transcript_id):
    '''
    transcript_id is index in df_transcript2
    translate CDSplus based on CDSplus and frame
    '''
    r = df_transcript2.loc[transcript_id]
    CDSplus = r['CDSplus']
    strand = r['strand']
    frame = int(r['frame'])
    CDSplus = Seq(CDSplus)
    if strand == '-':
        CDSplus = CDSplus.reverse_complement()
    CDSplus = CDSplus[frame:]
    CDSplus = CDSplus[:(len(CDSplus) // 3) * 3]
    AA_seq = str(CDSplus.translate().split('*')[0])
    if frame != 0:
        AA_seq = 'X' + AA_seq
    return AA_seq

# get loc for each nt in CDSplus
def getPostionsFromLocs(locs):
    '''
    locs like [(65564, 65573), (69036, 71585)]
    return a list of positions
    '''
    results = []
    for s,e in locs:
        results = results + list(range(s+1, e+1))#Note, there is a -1 operation for genomicLocs above
    return results

f_rc = lambda x:str(Seq(x).reverse_complement()) # get reverse completment of a sequence
def getMut(mutations, strand):
    '''
    given index of mutations in df_mutations,
    return a dataframe, with chr, pos, and mutations
    for substitution or insertion, no change for forward, and reverse completement for reverse strand
    for deletion, change so that ref include only one AA, alt change to empty
    '''
    tdf_mut = df_mutations.loc[mutations]
    tdf_mut['variant_id'] = tdf_mut.apply(lambda x:'{}-{}-{}-{}'.format(x['chr'], x['pos'], x['ref'], x['alt']),axis=1)
    
    results = []
    for row, r in tdf_mut.iterrows():
        chromosome = r['chr']
        pos = r['pos']
        ref = r['ref']
        alt = r['alt']
        variant_id = r['variant_id']
        if len(ref) == 1:
            if strand == '+':
                results.append([chromosome, pos, ref, alt, variant_id])
            else:
                results.append([chromosome, pos, f_rc(ref), f_rc(alt), variant_id])
        else:
            for n in range(1, len(ref)):
                p = pos + n
                r = ref[n]
                if strand == '+':
                    results.append([chromosome, p, r, '', variant_id])
                else:
                    results.append([chromosome, p, f_rc(r), '', variant_id])
    
    tdf_m = pd.DataFrame(results)
    tdf_m.columns = ['chr','pos','ref','alt','variant_id']
    return tdf_m

def getCodons(ttdf, AA_len=None):
    '''
    ttdf is with 'locs', 'bases', 'chr', 'strand'
    return a dataframe, in the form of (chromosome, strand, codon1_position, codon), like ('1','+',1512269, 'ATG')
    AA_len is the number of codons to keep. If None, return as long as possible. It is supporsed that the reading frame is correct and each line with only one base
    '''
    if AA_len is None:
        AA_len = ttdf.shape[0] // 3
    ttdf = ttdf.iloc[:AA_len * 3]
    tdf_result = ttdf.groupby(np.arange(AA_len * 3) // 3).agg({'chr':'first', 'strand':'first', 'locs':'first', 'bases':lambda x:''.join(x), 'variant_id':lambda x:','.join(list(dict.fromkeys([e for e in list(x) if pd.notnull(e)])))})
    tdf_result.columns = ['chr', 'strand', 'codon1', 'codon','variants']
    return tdf_result

def translateCDSplusWithMut(transcript_id):
    '''
    transcript_id is index in df_transcript3
    translate CDSplus based on CDSplus and frame, and mutations
    '''
#    t0 = time.time()
    r = df_transcript3.loc[transcript_id]
    locs = r['genomicLocs']
    locs = getPostionsFromLocs(locs)
    CDSplus = Seq(r['CDSplus'])
    AA_seq = r['AA_seq']
    AA_ori = AA_seq
    frame = int(r['frame'])
    strand = r['strand']
    chromosome = r['seqname']
    mutations = r['mutations']
    if strand == '-':
        CDSplus = CDSplus.reverse_complement()
        locs.reverse()
    if frame != 0:
        CDSplus = CDSplus[frame:]
        locs = locs[frame:]
        AA_seq = AA_seq[1:]
    df_CDSplus = pd.DataFrame()
    df_CDSplus['locs'] = locs
    df_CDSplus['bases'] = list(CDSplus)
    df_CDSplus['chr'] = chromosome
    df_CDSplus['strand'] = strand
#    print(time.time() - t0)
    tdf_m = getMut(mutations, strand)
#    print(time.time() - t0)
    # include mutation data
    df_CDSplus = df_CDSplus.merge(tdf_m, how='left', left_on = ['chr','locs', 'bases'], right_on=['chr','pos','ref'])
    # number of variant in CDSplus
#    n_variant_CDSplus = df_CDSplus[df_CDSplus['alt'].notnull()].shape[0]
    df_CDSplus['new_nt'] = df_CDSplus.apply(lambda x:x['alt'] if not pd.isnull(x['alt']) else x['bases'], axis=1)
#    print(time.time() - t0)
    # in two cases for gencode, the deletion is in the start codon, making the code not working. So if there is deletion in the start codon, making the first three new_nt unchanged
    if df_CDSplus.iloc[2]['new_nt'] != '':
        if df_CDSplus.iloc[0]['new_nt'] == '' or df_CDSplus.iloc[1]['new_nt'] == '':
            print(transcript_id, 'special case: start codon frame shift deletion')
        if df_CDSplus.iloc[0]['new_nt'] == '':
            df_CDSplus.loc[0, 'new_nt'] = df_CDSplus.iloc[0]['ref']
        if df_CDSplus.iloc[1]['new_nt'] == '':
            df_CDSplus.loc[1,'new_nt'] = df_CDSplus.iloc[1]['ref']
            
#    print(time.time() - t0)
    df_CDSref = df_CDSplus[['locs', 'bases', 'chr', 'strand','variant_id']].copy()
    df_CDSalt = df_CDSplus[['locs', 'new_nt', 'chr', 'strand','variant_id']].copy()
    # modify df_CDSalt so that each line with only one base
    df_CDSalt['bases'] = df_CDSalt['new_nt'].apply(list)
    df_CDSalt = df_CDSalt.explode('bases')
    df_CDSalt = df_CDSalt[df_CDSalt['bases'].notnull()]
#    print(time.time() - t0)
    # get codons from df_CDSref. codon is (chromosome, strand, start, codon), like ('1','+',1512269, 'ATG'). value is from AA_seq
    AA_len = len(AA_seq)
    AA_translate = str(CDSplus[:3*(len(CDSplus)//3)].translate().split('*')[0])
    if len(AA_translate) > AA_len:#in some cases, the protein sequences are not stop at stop codon but a little earlier
        AA_seq = AA_translate#if AA_translate shorter than AA_seq, there are usually non-standard codons, stop codon to AA
        AA_len = len(AA_seq)
    
    
    codons_ref = getCodons(df_CDSref, AA_len=AA_len)
    codons_alt = getCodons(df_CDSalt, AA_len=None)
    codons_ref['AA_ref'] = list(AA_seq)
    codons_ref['AA_index'] = list(range(1, AA_len+1))
    if frame !=0: codons_ref['AA_index'] = codons_ref['AA_index'] + 1
    # add ref AA to codon_alt
    codons_alt = codons_alt.merge(codons_ref, left_on=['chr','strand','codon1'], right_on=['chr','strand','codon1'], how='left')
    codons_alt.columns = ['chr','strand','codon1', 'codon_alt','variants', 'codon_ref', 'variants2','AA_ref','AA_index']
    codons_alt['AA_alt'] = codons_alt.apply(lambda x:x['AA_ref'] if x['codon_alt'] == x['codon_ref'] else str(Seq(x['codon_alt']).translate()), axis=1)
    codons_alt['AA_index'] = codons_alt['AA_index'].fillna(method='ffill')
    new_AA = str(''.join(codons_alt['AA_alt']).split('*')[0])
#    print(time.time() - t0)
    
    # check if final reading frame is the same
    t_last_frame = codons_alt[codons_alt['codon1'] == codons_ref.iloc[-1]['codon1']]
    if t_last_frame.shape[0] == 1:#final reading frame the same
        t_alt = codons_alt.loc[0:t_last_frame.index[0]]# cut to the last reading frame
        if '*' in set(t_alt['AA_alt']):
            t_alt = t_alt.loc[0:t_alt[t_alt['AA_alt'] == '*'].first_valid_index()] # cut to the first stop codon
            if pd.isnull(t_alt.iloc[-1]['AA_ref']):
                frameChange = True
            else:
                frameChange =False
        else:
            frameChange = False
    else:
        frameChange = True
    
#    print(time.time() - t0)
    # check if there is a stop gain
    t_alt = codons_alt.iloc[:len(AA_seq)]
    t_alt = t_alt[(t_alt['AA_alt'] == '*')]
    if t_alt.shape[0] > 0:
        stopGain = True
        x = t_alt.iloc[0].copy()
        if pd.isnull(x['AA_ref']):
            x['AA_ref'] = '-'
        AA_stopGain = x['AA_ref'] + str(int(x['AA_index'])) + x['AA_alt'] + "({})".format(x['variants'])
        stopAA_index = t_alt.iloc[0]['AA_index']
    else:
        stopGain = False
        AA_stopGain = ''
#    print(time.time() - t0)
    
    # check if there is a stop loss
    t_alt = codons_alt[codons_alt['codon1'] == codons_ref.iloc[-1]['codon1']]
    stopLoss = False
    stopLoss_pos = ''
    if t_alt.shape[0] > 0:
        if t_alt.index[0]+1 in codons_alt.index:
            if codons_alt.loc[t_alt.index[0]+1]['AA_alt'] != '*':
                stopLoss = True
                stopLoss_position = int(codons_alt.loc[t_alt.index[0]+1]['AA_index'])
                stopLoss_pos = str(stopLoss_position) + "({})".format(codons_alt.loc[t_alt.index[0]+1]['variants'])
                codons_alt = codons_alt.loc[:int(stopLoss_position-1)]
#    print(time.time() - t0)
    
    if stopGain:
        codons_alt = codons_alt.iloc[:len(new_AA)]#cut at valid sequences in new_AA
        codons_ref = codons_ref[codons_ref['AA_index'] < stopAA_index]
    
    if frameChange:
        codons_alt = codons_alt.loc[:codons_alt['AA_ref'].last_valid_index()]#cut at last match with codon_ref
        codons_ref = codons_ref[codons_ref['AA_index'] <= codons_alt['AA_index'].max()]#cut at no frame change codon with codon_alt
        
    codons_alt = codons_alt.iloc[:len(new_AA)]#cut at valid sequences in new_AA
        
    # count number of mutated AAs
    t_alt = codons_alt.dropna(subset=['AA_ref','AA_alt'])
    n_variant_AA = (t_alt['AA_ref'] != t_alt['AA_alt']).sum()
    if n_variant_AA > 0:
        t_alt = t_alt[t_alt['AA_ref'] != t_alt['AA_alt']]
        variant_AA = ';'.join(t_alt.apply(lambda x:x['AA_ref'] + str(int(x['AA_index'])) + x['AA_alt'] + '({})'.format(x['variants']), axis=1))
    else:
        variant_AA = ''
    # get deletions
    deletion_AA = set(codons_ref['AA_index']) - set(codons_alt['AA_index'])
    n_deletion_AA = len(deletion_AA)
    if n_deletion_AA > 0:
        t_ref = codons_ref[codons_ref['AA_index'].isin(deletion_AA)]
        deletion_AA = ';'.join(t_ref.apply(lambda x:x['AA_ref'] + str(int(x['AA_index'])) + '-' + '({})'.format(x['variants']), axis=1))
    else:
        deletion_AA = ''
    # get insertions
    t_alt = codons_alt[codons_alt['AA_ref'].isnull()]
    n_insertion_AA = t_alt.shape[0]
    if n_insertion_AA > 0:
        insertion_AA = ';'.join(t_alt.apply(lambda x:'-' + str(int(x['AA_index'])) + x['AA_alt'] + '({})'.format(x['variants']), axis=1))
    else:
        insertion_AA = ''
        
    if frame != 0:
        new_AA = 'X' + new_AA
#    print(time.time() - t0)
    if not frameChange:
        if n_variant_AA + n_deletion_AA + n_insertion_AA == 0 and (not stopGain):
            if AA_ori in new_AA:#new_AA is longer than AA_ori at the end of the sequence, stop loss
                new_AA = new_AA #AA_ori
            elif new_AA in AA_translate:
                print(transcript_id, 'AA_translate longer than provided AA sequence, the gtf file not consistent with the protein file')
            else:
                print(transcript_id,'something wrong?')
                print(AA_ori)
                print(new_AA)
    
    return new_AA, frameChange, stopGain, AA_stopGain, stopLoss, stopLoss_pos, n_variant_AA, n_deletion_AA, n_insertion_AA, variant_AA, insertion_AA, deletion_AA

def translateCDSplusWithMut2(transcript_id):
    try:
        return translateCDSplusWithMut(transcript_id)
    except:
        print(transcript_id, 'cannot be processed properly, please check')
        return transcript_id


description = '''output a new reference protein set by with the variants data
file_genome: the reference genome sequence in fasta format
file_proteins: protein sequences in fasta format. Can be None. If None, translate the protein sequences based on the CDS annotation in the gtf file
file_mutations: a file stores the variants
outprefix: the prefix of output files
cpu_counts: number of CPUs to use when running the program
'''

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-g','--genome', help = 'the reference genome sequence in fasta format. It can be a gzip file', required=True)
    parser.add_argument('-p','--protein', help = 'protein sequences in fasta format. It can be a gzip file. default None. If None, translate the protein sequences based on the CDS annotation in the gtf file', default=None)
    parser.add_argument('-f', '--gtf', help='gtf file with CDS and exon annotations', required=True)
    parser.add_argument('-m', '--mutations', help='a file stores the variants', required=True)
    parser.add_argument('-t', '--threads', help='number of threads/CPUs to run the program. default, use all CPUs available', type=int, default=os.cpu_count())
    parser.add_argument('-o', '--out', help='output prefix. two file will be output. One is the annotation for mutated transcripts, one is the protein sequences. {out}_mutations.csv, {out}_protein.fa')
    
    f = parser.parse_args()
    
    file_genome = f.genome
    file_proteins = f.protein
    file_gtf = f.gtf
    file_mutations = f.mutations
    outprefix = f.out
    cpu_counts = os.cpu_count()
    
    
    # parse mutation file
    df_mutations = parse_mutations(file_mutations)
    df_transcript2 = get_df_transcript2(file_gtf, file_proteins, file_genome, cpu_counts)
    
    # assign mutations to each transcript. 
    # add column 'mutations' to df_transcript2, which stores index on the mutations in df_mutations
    # add column 'mutations', with list of mutation index in df_mutations
    pool = Pool(cpu_counts)
    results = pool.map(getMutations, df_transcript2.index)
    pool.close()
    pool.join()
    df_transcript2['mutations'] = results
    
    # add column 'AA_translate'
    pool = Pool(cpu_counts)
    df_transcript2['AA_translate'] = pool.map(translateCDSplus, df_transcript2.index)
    pool.close()
    pool.join()
    
    # split df_transcript2 to with mutations and no_mutations
    df_transcript_noMut = df_transcript2[df_transcript2['mutations'].apply(lambda x:len(x)==0)]#57677
    df_transcript3 = df_transcript2[df_transcript2['mutations'].apply(lambda x:len(x)!=0)].copy()#53132
    
    # save df_transcript_noMut
    ## if with AA_seq, save directly
    df_transcript_noMutWithAA = df_transcript_noMut[df_transcript_noMut['AA_seq'].notnull()]#47459
    df_transcript_noMutNoAA = df_transcript_noMut[df_transcript_noMut['AA_seq'].isnull()]#10218
    
    print('transcripts with no mutations:', df_transcript_noMut.shape[0])
    print('transcripts with no mutations, with protein sequence:', df_transcript_noMutWithAA.shape[0],', AA will be unchanged')
    print('transcripts with no mutations, without protein sequence:', df_transcript_noMutNoAA.shape[0], ', AA is translated from the CDS frame till the first stop codon')
    
    # add AA_seq if with no AA_seq for transcript_ids without AA_seq
    print('transcripts with mutations in CDSplus:', df_transcript3.shape[0])
    print('transcripts with mutations, {n} without protein sequence. translated from CDSplus seq'.format(n=df_transcript3.AA_seq.isnull().sum()))
    # add AA_seq for df_transcript3 if with no AA_seq
    df_transcript3['AA_seq'] = df_transcript3.apply(lambda x:x['AA_seq'] if pd.notnull(x['AA_seq']) else x['AA_translate'], axis=1)
    
    pool = Pool(cpu_counts)
    results = pool.map(translateCDSplusWithMut2, list(df_transcript3.index))
    pool.close()
    pool.join()
    tdf = pd.DataFrame(results)
    tdf.columns = 'new_AA frameChange  stopGain  AA_stopGain  stopLoss  stopLoss_pos  n_variant_AA  n_deletion_AA  n_insertion_AA  variant_AA  insertion_AA  deletion_AA'.split()
    for col in tdf.columns:
        df_transcript3[col] = list(tdf[col])
    
    # get sequence length of ref_AA and alt_AA
    df_transcript3['len_ref_AA'] = df_transcript3['AA_seq'].str.len()
    df_transcript3['len_alt_AA'] = df_transcript3['new_AA'].str.len()
    
    # save mutation annotation
    df_sum_mutations = df_transcript3[['seqname','strand','gene_id','protein_id','hgnc_id','gene_name','frameChange','stopGain', 'AA_stopGain', 'stopLoss', 'stopLoss_pos', 'n_variant_AA', 'n_deletion_AA', 'n_insertion_AA', 'variant_AA', 'insertion_AA', 'deletion_AA', 'len_ref_AA', 'len_alt_AA']]
    # only keep transcripts with mutations
    df_sum_mutations = df_sum_mutations[df_sum_mutations['frameChange'] | df_sum_mutations['stopGain'] | df_sum_mutations['stopLoss'] | (df_sum_mutations['n_variant_AA'] > 0) | (df_sum_mutations['n_deletion_AA'] > 0) | (df_sum_mutations['n_insertion_AA'] > 0)]
    
    print('of {n1} transcripts from {n2} genes with mutations in CDSplus regions, {n3} transcripts from {n4} genes with protein sequences altered finally'.format(n1 = df_transcript3.shape[0], n2 = df_transcript3['gene_id'].unique().shape[0], n3 = df_sum_mutations.shape[0], n4 = df_sum_mutations['gene_id'].unique().shape[0]))
    
    # print some summary
    print('total frame change:', df_sum_mutations.frameChange.sum())
    print('total stop gain:', df_sum_mutations.stopGain.sum())
    print('total stop loss:', df_sum_mutations.stopLoss.sum())
    print('total AA mutations:', df_sum_mutations.n_variant_AA.sum())
    print('total AA deletions:', df_sum_mutations.n_deletion_AA.sum())
    print('total AA insertions:', df_sum_mutations.n_insertion_AA.sum())
    
    outfilename = outprefix +'_mutations.csv'
    if not os.path.exists(os.path.dirname(outfilename)):
        os.makedirs(os.path.dirname(outfilename))
    df_sum_mutations.to_csv(outfilename, sep='\t')
    
    # save protein sequences
    ## if with AA_seq, save directly
    df_save1 = df_transcript_noMutWithAA.copy()#47459
    df_save2 = df_transcript_noMutNoAA.copy()#10218
    df_save3 = df_transcript3[~(df_transcript3['frameChange'] | df_transcript3['stopGain'] | df_transcript3['stopLoss'] | (df_transcript3['n_variant_AA'] > 0) | (df_transcript3['n_deletion_AA'] > 0) | (df_transcript3['n_insertion_AA'] > 0))].copy()
    df_save4 = df_transcript3.loc[df_sum_mutations.index].copy()
    df_save1.reset_index(inplace=True)
    df_save2.reset_index(inplace=True)
    df_save3.reset_index(inplace=True)
    df_save4.reset_index(inplace=True)
    
    fout = open(outprefix + '_protein.fa','w')
    if df_save1.shape[0] > 0:
        fout.write(''.join(df_save1.apply(lambda x:'>{}|{}|{}|{}|{} noMut_withAA\n{}\n'.format(x['transcript_id'], x['protein_id'], x['gene_id'], x['hgnc_id'], x['gene_name'],x['AA_seq']), axis=1)))
    if df_save2.shape[0] > 0:
        fout.write(''.join(df_save2.apply(lambda x:'>{}|{}|{}|{}|{} noMut_withoutAA\n{}\n'.format(x['transcript_id'], x['protein_id'], x['gene_id'], x['hgnc_id'], x['gene_name'],x['AA_translate']), axis=1)))
    if df_save3.shape[0] > 0:
        fout.write(''.join(df_save3.apply(lambda x:'>{}|{}|{}|{}|{} withMut_noChange\n{}\n'.format(x['transcript_id'], x['protein_id'], x['gene_id'], x['hgnc_id'], x['gene_name'],x['AA_seq']), axis=1)))
    if df_save4.shape[0] > 0:
        fout.write(''.join(df_save4.apply(lambda x:'>{}|{}|{}|{}|{} withMut_AAchange\n{}\n'.format(x['transcript_id'], x['protein_id'], x['gene_id'], x['hgnc_id'], x['gene_name'],x['new_AA']), axis=1)))
    
    fout.close()