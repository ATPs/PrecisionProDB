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


def openFile(filename):
    '''open text or gzipped text file
    '''
    if filename.endswith('.gz'):
        return gzip.open(filename,'rt')
    return open(filename)


def parse_proteins(file_proteins):
    '''
    return a dataframe with file_proteins. file_proteins is processed by perGeno, with '\t' to separate protein_id and original description
    '''
    # get protein sequences
    if file_proteins.endswith('.gz'):
        tls =list(SeqIO.parse(gzip.open(file_proteins,'rt'), 'fasta'))
    else:
        tls =list(SeqIO.parse(open(file_proteins,'r'), 'fasta'))
    df_protein = pd.DataFrame([[seq.id, seq.description.split('\t',maxsplit=1)[1], str(seq.seq).strip('*')] for seq in tls])#remvoe '*' which represent stop codons
    df_protein.columns = ['protein_id','protein_description', 'AA_seq']
    print('from the protein file, totally', len(tls), 'protein sequences.')
    return df_protein

def parse_genome(file_genome):
    '''
    file_genome can be a gzip file. It should be a fasta file
    given a genome file, return a string of sequences
    '''
    chromsome_seq = SeqIO.read(openFile(file_genome), 'fasta')
    return str(chromsome_seq.seq).upper()


def getCDSplus(tdf):
    '''
    tdf is (transcript_id, tdf)
    given a transcript_id, get the DNA_seq based on the gtf annotation
    dcgtf_transcript is a dict with transcript_id as key, and gtf dataframe for that transcript as value
    '''
    transcript_id, tdf = tdf
    tdf = tdf[tdf['feature'].isin(['exon','CDS'])]
    strand = tdf['strand'].iloc[0]
    if strand =='+':
        tdf = tdf.sort_values(by=['start','end'])# sort by location from small to large if positive strand
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
    # cases that only with CDS annotations
    if 'exon' not in set(tdf['feature']):
        tlocs = [[i-1,j] for i,j in zip(tdf['start'], tdf['end'])]# all locs of CDS
    else:
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
        tlocs = [[i-1,j] for i,j in zip(tdf['start'], tdf['end'])]

    tlocs = sorted(tlocs, key=lambda x:x[0]) # sort by location from small to large
    return transcript_id, frame, tlocs
    

def getCDSstrand(CDS, AA_seq):
    '''
    return CDS strand. if CDS translation equal AA_seq, return '+'
    if reverse complement of CDS translated to AA_seq, return '-'
    otherwise, print warning, return '.'
    '''
    CDS_f = Seq(CDS)
    CDS_r = CDS_f.reverse_complement()
    AA_f = str(CDS_f.translate().split('*')[0])
    AA_r = str(CDS_r.translate().split('*')[0])
    if AA_f == AA_seq:
        return '+'
    if AA_r == AA_seq:
        return '-'
    
    return '.'

def getStrand(r):
    '''
    r is a row in df_transcript2
    return the correct strand based on CDS sequence and AA_seq
    '''
    AA_seq = r['AA_seq']
    CDS = r['CDS']
    strand = r['strand']
    strand_n = getCDSstrand(CDS, AA_seq)
    if strand_n == '.':
        print(r.name,'something wrong with CDS region')
    elif strand_n != strand:
        print(r.name, 'strand changed from {} to {}'.format(strand, strand_n))
    else:
        pass
    return strand_n

def get_df_transcript2(file_gtf, file_proteins, file_genome, cpu_counts):
    '''
    get df_transcripts, add AA_seq, add CDSplus and other notes
    for some reason, the strand of sequences in the gtf file is not accurate. need to calculate
    '''
    #global datatype
    # get df_gtf
    df_gtf = pd.read_csv(file_gtf, sep='\t',header=None, comment='#')
    df_gtf.columns = ['seqname','source','feature','start','end','score','strand','frame','protein_id']
    # change df_gtf seqname to str
    df_gtf['seqname'] = df_gtf['seqname'].astype(str)
    
    # get df_protein
    df_protein = parse_proteins(file_proteins)
    
    # get genome DNA sequences
    chromosome_seq = parse_genome(file_genome)
    print('finish reading genome file')
    
    # get protein_id in df_gtf
    protein_ids_gtf = set(df_gtf[df_gtf['protein_id'].notnull()]['protein_id'])
    print('number of proteins from the gtf file', len(protein_ids_gtf)) 
    
    # some protein_id only found in the protein fasta file
    df_transcript2 = df_protein.set_index('protein_id').copy()
    # add strand
    dc_gtfpr2strand = dict(zip(df_gtf['protein_id'], df_gtf['strand']))
    df_transcript2['strand'] = [dc_gtfpr2strand[e] for e in df_transcript2.index]

    df_cds = df_gtf[df_gtf['feature'] == 'CDS']
    tdc_pr2cds = {k:v for k,v in df_cds.groupby('protein_id')}
    
    # add CDSLoc and CDS sequence
    dc_pr2cdsLocs = {}
    for k,v in tdc_pr2cds.items():
        v = v.sort_values(by='start')
        cdsLocs = [[i-1,j] for i,j in zip(v['start'],v['end'])]
        dc_pr2cdsLocs[k] = cdsLocs
    
    df_transcript2['CDSloc'] = [dc_pr2cdsLocs[e] for e in df_transcript2.index]
    df_transcript2['CDS'] = df_transcript2.apply(lambda x:''.join([chromosome_seq[i:j] for i,j in x['CDSloc']]), axis=1)
    
    # correct strand for gtf format data from transdecoder
    if datatype == 'gtf':
        strand_n = df_transcript2.apply(getStrand, axis=1)
        df_transcript2['strand'] = strand_n

    # get CDSplus sequences
    df_gtf_group = {k:v for k,v in df_gtf.groupby('protein_id') }
    
    pool = Pool(cpu_counts)
    results = pool.map(getCDSplus, df_gtf_group.items())
    pool.close()
    pool.join()
    
    
    tdf = pd.DataFrame(results)
    tdf.columns = ['protein_id','frame', 'genomicLocs']
    df_transcript2 = df_transcript2.merge(tdf, left_index=True, right_on='protein_id', how='left')
    print('finishing get locs and frame')
    df_transcript2['genomicStart'] = df_transcript2['genomicLocs'].apply(lambda x:None if not isinstance(x,list) else x[0][0])
    df_transcript2['genomicEnd'] = df_transcript2['genomicLocs'].apply(lambda x:None if not isinstance(x,list) else x[-1][1])
    
    df_transcript2['CDSplus'] = df_transcript2.apply(lambda x:''.join([chromosome_seq[i:j] for i,j in x['genomicLocs']]), axis=1)
    
    df_transcript2 = df_transcript2.set_index('protein_id')
    return df_transcript2

def getMutations(transcript_id):
    '''
    given a row in df_transcript2 of transcript_id, return the index of mutations where it locates inside genomicStart and genomicEnd
    Note the genomicStart should be add 1, as it is index starting from 0
    '''
    r = df_transcript2.loc[transcript_id]
    genomicLocs = r['genomicLocs']
    results = []
    tdf = df_mutations
    for genomicStart, genomicEnd in genomicLocs:
        genomicStart = genomicStart + 1
        tdf = tdf[(tdf['pos'] >= genomicStart) & (tdf['pos'] <= genomicEnd) & (tdf['pos_end'] >= genomicStart) & (tdf['pos_end'] <= genomicEnd)]
        results += list(tdf.index)
    return results

def getMutations_faster(transcript_id):
    '''
    given a row in df_transcript2 of transcript_id, return the index of mutations where it locates inside genomicStart and genomicEnd
    Note the genomicStart should be add 1, as it is index starting from 0
    '''
    r = df_transcript2.loc[transcript_id]
    genomicLocs = r['genomicLocs']
    results = []
    positions = df_mutations['pos']
    for genomicStart, genomicEnd in genomicLocs:
        genomicStart = genomicStart + 1
        insert_start = positions.searchsorted(genomicStart,side='left')
        insert_end = positions.searchsorted(genomicEnd, side='right')
        positions_keep = positions.iloc[insert_start:insert_end]
        results += list(positions_keep.index)
    
    # also consider the end of the mutation
    positions2 = df_mutations['pos_end']
    results2 = []
    for genomicStart, genomicEnd in genomicLocs:
        genomicStart = genomicStart + 1
        insert_start = positions2.searchsorted(genomicStart,side='left')
        insert_end = positions2.searchsorted(genomicEnd, side='right')
        positions_keep = positions2.iloc[insert_start:insert_end]
        results2 += list(positions_keep.index)
    
    return [e for e in results if e in results2]

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
    #global chromosome
    tdc_result = {}
    r = df_transcript3.loc[transcript_id]
    locs = r['genomicLocs']
    locs = getPostionsFromLocs(locs)
    CDSplus = Seq(r['CDSplus'])
    AA_seq = r['AA_seq']
    AA_ori = AA_seq
    AA_translate = r['AA_translate']
    frame = int(r['frame'])
    strand = r['strand']
    mutations = r['mutations']
    tdf_m = getMut(mutations, strand)
    
    if strand == '-':
        CDSplus = CDSplus.reverse_complement()
        locs.reverse()
    if frame != 0:
        CDSplus = CDSplus[frame:]
        locs = locs[frame:]
        AA_seq = AA_seq[1:]
        AA_translate = AA_translate[1:]

    AA_len = len(AA_seq)

    # check if with non-standard stop codon. keep the non-standard stop codon by trim CDSplus and locs
    if AA_len < len(AA_translate) and len(CDSplus)/3 - AA_len > 1:#translation stop not at stop codon
        tdc_result['nonStandardStopCodon'] = 1
        CDSplus = CDSplus[:3*AA_len]
        locs = locs[:3*AA_len]
    
    df_CDSplus = pd.DataFrame()
    df_CDSplus['locs'] = locs
    df_CDSplus['bases'] = list(CDSplus)
    df_CDSplus['chr'] = chromosome
    df_CDSplus['strand'] = strand
    # include mutation data
    df_CDSplus = df_CDSplus.merge(tdf_m, how='left', left_on = ['chr','locs', 'bases'], right_on=['chr','pos','ref'])
    # number of variant in CDSplus
    df_CDSplus['new_nt'] = df_CDSplus.apply(lambda x:x['alt'] if not pd.isnull(x['alt']) else x['bases'], axis=1)
    # in two cases for gencode, the deletion is in the start codon, making the code not working. So if there is deletion in the start codon, making the first three new_nt unchanged
    if df_CDSplus.iloc[2]['new_nt'] != '':
        if df_CDSplus.iloc[0]['new_nt'] == '' or df_CDSplus.iloc[1]['new_nt'] == '':
            print(transcript_id, 'special case: start codon frame shift deletion')
        if df_CDSplus.iloc[0]['new_nt'] == '':
            df_CDSplus.loc[0, 'new_nt'] = df_CDSplus.iloc[0]['ref']
        if df_CDSplus.iloc[1]['new_nt'] == '':
            df_CDSplus.loc[1,'new_nt'] = df_CDSplus.iloc[1]['ref']
    
    # if no special condon, translate directly
    if AA_seq in AA_translate:
        CDSnew = ''.join(df_CDSplus['new_nt'])
        new_AA = str(Seq(CDSnew[:3*(len(CDSnew) // 3)]).translate().split('*')[0])
        tdc_result['new_AA'] = new_AA
        # no AA changed, return original sequence
        if new_AA == AA_seq:
            if frame != 0:
                tdc_result['new_AA'] = 'X' + tdc_result['new_AA']
            return tdc_result
        # only AA change
        if len(new_AA) == len(AA_seq):
            t_alt = pd.DataFrame()
            for n in range(len(AA_seq)):
                if AA_seq[n] != new_AA[n]:
                    t_alt.loc[n, 'AA_ref'] = AA_seq[n]
                    t_alt.loc[n, 'AA_alt'] = new_AA[n]
                    t_alt.loc[n, 'AA_index'] = n + 1
                    t_alt.loc[n, 'variants'] = ','.join(df_CDSplus.loc[n*3:(n+1)*3,'variant_id'].dropna().drop_duplicates())
            tdc_result['n_variant_AA'] = t_alt.shape[0]
            tdc_result['variant_AA'] = ';'.join(t_alt.apply(lambda x:x['AA_ref'] + str(int(x['AA_index'])) + x['AA_alt'] + '({})'.format(x['variants']), axis=1))
            if frame != 0:
                tdc_result['new_AA'] = 'X' + tdc_result['new_AA']
            return tdc_result
            
    # translation with non-standard codons or complicate cases
    df_CDSref = df_CDSplus[['locs', 'bases', 'chr', 'strand','variant_id']].copy()
    df_CDSalt = df_CDSplus[['locs', 'new_nt', 'chr', 'strand','variant_id']].copy()
    # modify df_CDSalt so that each line with only one base
    df_CDSalt['bases'] = df_CDSalt['new_nt'].apply(list)
    df_CDSalt = df_CDSalt.explode('bases')
    df_CDSalt = df_CDSalt[df_CDSalt['bases'].notnull()]
    # get codons from df_CDSref. codon is (chromosome, strand, start, codon), like ('1','+',1512269, 'ATG'). value is from AA_seq
    
    codons_ref = getCodons(df_CDSref, AA_len=AA_len)
    codons_alt = getCodons(df_CDSalt, AA_len=None)
    codons_ref['AA_ref'] = list(AA_seq)
    codons_ref['AA_index'] = list(range(1, AA_len+1))
    if frame !=0: codons_ref['AA_index'] = codons_ref['AA_index'] + 1
    # add ref AA to codon_alt
    codons_alt = codons_alt.merge(codons_ref, left_on=['chr','strand','codon1'], right_on=['chr','strand','codon1'], how='left')
    codons_alt.columns = ['chr','strand','codon1', 'codon_alt','variants', 'codon_ref', 'variants_ref','AA_ref','AA_index']
    codons_alt.loc[codons_alt.duplicated('AA_index'),['codon_ref','AA_ref','AA_index']] = np.NaN # only use each AA_ref Once
    codons_alt['AA_alt'] = codons_alt.apply(lambda x:x['AA_ref'] if x['codon_alt'] == x['codon_ref'] else str(Seq(x['codon_alt']).translate()), axis=1)
    codons_alt['AA_index'] = codons_alt['AA_index'].fillna(method='ffill')
    new_AA = str(''.join(codons_alt['AA_alt']).split('*')[0])
    
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
    
    # check if there is a stop gain
    t_alt = codons_alt.iloc[:len(AA_seq)]
    t_alt = t_alt[(t_alt['AA_alt'] == '*') & ((t_alt['variants'] != '') | t_alt['variants_ref'].notnull())]
    if t_alt.shape[0] > 0:
        stopGain = True
        x = t_alt.iloc[0].copy()
        if pd.isnull(x['AA_ref']):
            x['AA_ref'] = '-'
        if x['variants'] == '' and pd.notnull(x['variants_ref']):
            x['variants'] = x['variants_ref']
        AA_stopGain = x['AA_ref'] + str(int(x['AA_index'])) + x['AA_alt'] + "({})".format(x['variants'])
        stopAA_index = t_alt.iloc[0]['AA_index']
    else:
        stopGain = False
        AA_stopGain = ''
    
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
    
    if stopGain:
        codons_alt = codons_alt.iloc[:len(new_AA)]#cut at valid sequences in new_AA
        codons_ref = codons_ref[codons_ref['AA_index'] < stopAA_index]
    
    if frameChange:
        codons_alt = codons_alt.loc[:codons_alt['AA_ref'].last_valid_index()]#cut at last match with codon_ref
        codons_ref = codons_ref[codons_ref['AA_index'] <= codons_alt['AA_index'].max()]#cut at no frame change codon with codon_alt
        
    codons_alt = codons_alt.iloc[:len(new_AA)]#cut at valid sequences in new_AA
        
    # count number of mutated AAs
    t_alt = codons_alt.dropna(subset=['AA_ref','AA_alt'])
    t_alt = t_alt[(t_alt['AA_ref'] != t_alt['AA_alt']) & (t_alt['variants'] != '')]
    n_variant_AA = t_alt.shape[0]
    if n_variant_AA > 0:
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
    
    tdc_result['new_AA'] = new_AA
    tdc_result['frameChange'] = frameChange
    tdc_result['stopGain'] = stopGain
    tdc_result['AA_stopGain'] = AA_stopGain
    tdc_result['stopLoss'] = stopLoss
    tdc_result['stopLoss_pos'] = stopLoss_pos
    tdc_result['n_variant_AA'] = n_variant_AA
    tdc_result['n_deletion_AA'] = n_deletion_AA
    tdc_result['n_insertion_AA'] = n_insertion_AA
    tdc_result['variant_AA'] = variant_AA
    tdc_result['insertion_AA'] = insertion_AA
    tdc_result['deletion_AA'] = deletion_AA
    return tdc_result

def translateCDSplusWithMut2(transcript_id):
    try:
        return translateCDSplusWithMut(transcript_id)
    except:
        print(transcript_id, 'cannot be processed properly, please check')
        return {}

def checkIfAAtranslatedFromGenome(r):
    '''given a r is a row in df_transcript2, return if the protein sequence agrees with the genome annotation
    '''
    if len(r['AA_seq']) == len(r['AA_translate']):
        return True

    CDS = r['CDS']
    strand = r['strand']
    frame = int(r['frame'])
    CDS = Seq(CDS)
    if strand == '-':
        CDS = CDS.reverse_complement()
    CDS = CDS[frame:]
    CDS = CDS[:(len(CDS) // 3) * 3]
    CDS_len = len(CDS)
    if str(CDS[-3:]) in ['TAA','TGA','TAG']:
        CDS_len -= 3
    if CDS_len / len(r['AA_seq']) == 3:
        return True
    return False




description = '''output a new reference protein set by with the variants data
'''

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-g','--genome', help = 'the reference genome sequence in fasta format. It can be a gzip file', required=True)
    parser.add_argument('-p','--protein', help = 'protein sequences in fasta format. It can be a gzip file.', default=None)
    parser.add_argument('-f', '--gtf', help='gtf file with CDS and exon annotations', required=True)
    parser.add_argument('-m', '--mutations', help='a file stores the variants', required=True)
    parser.add_argument('-t', '--threads', help='number of threads/CPUs to run the program. default, use all CPUs available', type=int, default=os.cpu_count())
    parser.add_argument('-o', '--out', help='output prefix. two file will be output. One is the annotation for mutated transcripts, one is the protein sequences. {out}.aa_mutations.csv, {out}.mutated_protein.fa')
    parser.add_argument('-c', '--chromosome', help='''chromosome name/id, default="chr1" ''', default='chr1', type=str)
    parser.add_argument('-d', '--datatype', help='''input datatype. choice of 'GENCODE_GTF','GENCODE_GFF3', 'RefSeq', 'gtf', default 'GENCODE_GTF' ''', type=str, choices=['GENCODE_GTF','GENCODE_GFF3', 'RefSeq', 'gtf'], default='GENCODE_GTF')
    
    f = parser.parse_args()
    
    file_genome = f.genome
    file_proteins = f.protein
    file_gtf = f.gtf
    file_mutations = f.mutations
    outprefix = f.out
    cpu_counts = f.threads
    chromosome = f.chromosome
    datatype = f.datatype
    
    
    # parse mutation file
    df_mutations = pd.read_csv(file_mutations, sep='\t',low_memory=False)
    df_mutations['pos_end'] = df_mutations['pos'] + df_mutations['ref'].str.len() -1 # add a column, 'pos_end' to include all span of the mutation in reference genome
    df_mutations = df_mutations.sort_values(by = 'pos')
    df_mutations['chr'] = chromosome
    
    df_transcript2 = get_df_transcript2(file_gtf, file_proteins, file_genome, cpu_counts)
    df_transcript2['seqname'] = chromosome
    
    # assign mutations to each transcript. 
    # add column 'mutations' to df_transcript2, which stores index on the mutations in df_mutations
    # add column 'mutations', with list of mutation index in df_mutations
    pool = Pool(cpu_counts)
    results = pool.map(getMutations_faster, df_transcript2.index)
    pool.close()
    pool.join()
    df_transcript2['mutations'] = results
    
    # add column 'AA_translate'
    pool = Pool(cpu_counts)
    df_transcript2['AA_translate'] = pool.map(translateCDSplus, df_transcript2.index)
    pool.close()
    pool.join()
    
    # for RefSeq, sometimes the GTF annotation does not agree with the protein sequence, which means the len(CDS)/len(protein) != 3, find those cases and do not change the proteins
    if datatype == 'RefSeq':
        tdf_special = df_transcript2.loc[~df_transcript2.apply(checkIfAAtranslatedFromGenome,axis=1)]
        print('chromosome', chromosome, list(tdf_special.index), 'not translated from the CDS sequences in the genome. do not change')
        for transcript_id in tdf_special.index:
            df_transcript2.at[transcript_id, 'mutations'] = []
    
    # split df_transcript2 to with mutations and no_mutations
    df_transcript_noMut = df_transcript2[df_transcript2['mutations'].apply(lambda x:len(x)==0)]
    df_transcript3 = df_transcript2[~df_transcript2.index.isin(set(df_transcript_noMut.index))].copy()
    
    # save df_transcript_noMut
    print('transcripts with no mutations:', df_transcript_noMut.shape[0])
    
    # add AA_seq if with no AA_seq for transcript_ids without AA_seq
    print('transcripts with mutations in CDSplus:', df_transcript3.shape[0])
    
    pool = Pool(cpu_counts)
    results = pool.map(translateCDSplusWithMut2, list(df_transcript3.index))
    pool.close()
    pool.join()
    tdf = pd.DataFrame(results)
    for col in tdf.columns:
        df_transcript3[col] = list(tdf[col])
    
    # get sequence length of ref_AA and alt_AA
    df_transcript3['len_ref_AA'] = df_transcript3['AA_seq'].str.len()
    if 'new_AA' in df_transcript3.columns:
        df_transcript3['len_alt_AA'] = df_transcript3['new_AA'].str.len()
    
    # save mutation annotation
    columns_keep = ['seqname','strand','frameChange','stopGain', 'AA_stopGain', 'stopLoss', 'stopLoss_pos', 'nonStandardStopCodon', 'n_variant_AA', 'n_deletion_AA', 'n_insertion_AA', 'variant_AA', 'insertion_AA', 'deletion_AA', 'len_ref_AA', 'len_alt_AA']
    columns_keep = [e for e in columns_keep if e in df_transcript3.columns]
    df_sum_mutations = df_transcript3[df_transcript3['AA_seq'] != df_transcript3['new_AA']][columns_keep]
    
    outfilename = outprefix +'.aa_mutations.csv'
    if not os.path.exists(os.path.dirname(outfilename)):
        os.makedirs(os.path.dirname(outfilename))
    df_sum_mutations.to_csv(outfilename, sep='\t')
    
    # save proteins
    fout = open(outprefix + '.mutated_protein.fa','w')
    for protein_id, r in df_transcript_noMut.iterrows():
        fout.write('>{}\n{}\n'.format(r['protein_description'], r['AA_seq']))
    for protein_id, r in df_transcript3.iterrows():
        if pd.notnull(r['new_AA']):
            fout.write('>{}\n{}\n'.format(r['protein_description'], r['new_AA']))
        else:
            fout.write('>{}\n{}\n'.format(r['protein_description'], r['AA_seq']))
    fout.close()