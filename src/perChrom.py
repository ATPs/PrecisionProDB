import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import gzip
import numpy as np
from multiprocessing import Pool
import os
import time
import numpy as np
import re

def openFile(filename):
    '''open text or gzipped text file
    '''
    if filename.endswith('.gz'):
        return gzip.open(filename,'rt')
    return open(filename)

def f_rc(x):
    '''get reverse completment of a sequence
    '''
    return str(Seq(x).reverse_complement())


def parse_proteins(file_proteins):
    '''
    return a dataframe with file_proteins. file_proteins is processed by perGeno, with '\t' to separate protein_id and original description
    '''
    # get protein sequences
    if file_proteins.endswith('.gz'):
        tls =list(SeqIO.parse(gzip.open(file_proteins,'rt'), 'fasta'))
    else:
        tls =list(SeqIO.parse(open(file_proteins,'r'), 'fasta'))
    df_protein = pd.DataFrame([[seq.id, seq.description.split('\t',maxsplit=1)[1], str(seq.seq).strip('*')] for seq in tls])#remove '*' which represent stop codons
    df_protein.columns = ['protein_id','protein_description', 'AA_seq']
    print('from the protein file, totally', len(tls), 'protein sequences.')
    df_protein['protein_id_fasta'] = df_protein['protein_description'].apply(lambda x:x.split()[0])
    df_protein['protein_anno'] = df_protein.apply(lambda x:x['protein_description'].split(maxsplit=1)[-1] if x['protein_description'].split(maxsplit=1)[-1] != x['protein_id_fasta'] else '', axis=1)
    return df_protein



def parse_mutation(file_mutations, chromosome=None):
    '''
    Parse mutation file.
    
    This function parses a mutation file, which can be either a string representing a single mutation or a file containing multiple mutations. 
    If the input is a string like "1-55051215-G-GA" or "chr1-55051215-G-GA" , it uses a regular expression to extract the chromosome, position, reference allele, and alternate allele. If the input is a file, it reads the file into a DataFrame and processes it accordingly. The function also handles the calculation of the end position for each mutation and optionally assigns a specific chromosome if provided.
    
    Parameters:
    file_mutations (str): The path to the mutation file or a string representing a single mutation.
    chromosome (str, optional): The chromosome to assign to the mutations. Defaults to None.
    
    Returns:
    pd.DataFrame: A DataFrame containing the parsed mutations with columns 'chr', 'pos', 'ref', 'alt', and 'pos_end'.
    '''
    # Use regex to parse the string
    pattern = re.compile(r'(chr)?(\d+)-(\d+)-([A-Za-z]+)-([A-Za-z]+)')
    match = pattern.match(file_mutations)
    if match and (not os.path.exists(file_mutations)):
        print('file_mutations is a single mutation', file_mutations)
        chr_name = match.group(2)
        pos = int(match.group(3))
        ref = match.group(4)
        alt = match.group(5)
        df_mutations = pd.DataFrame({
            'chr': [chr_name],
            'pos': [pos],
            'ref': [ref],
            'alt': [alt]
        })
        df_mutations['pos_end'] = df_mutations['pos'] + df_mutations['ref'].str.len() - 1
        if chromosome:
            df_mutations['chr'] = chromosome
    else:
        df_mutations = pd.read_csv(file_mutations, sep='\t', low_memory=False)
        df_mutations['pos_end'] = df_mutations['pos'] + df_mutations['ref'].str.len() - 1
        df_mutations = df_mutations.sort_values(by='pos')
        if chromosome:
            df_mutations['chr'] = chromosome
    
    df_mutations = df_mutations.sort_values(by=['chr','pos'])
    return df_mutations


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
    '''get CDS strand based on AA_seq
    Here CDS is always the forward strand
    if CDS translation equal AA_seq, return '+'
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


# get loc for each nt in CDSplus
def getPositionsFromLocs(locs):
    '''
    locs like [(65564, 65573), (69036, 71585)]
    return a list of positions
    '''
    results = []
    for s,e in locs:
        results = results + list(range(s+1, e+1))#Note, there is a -1 operation for genomicLocs above
    return results

def getPositionsFromLocs_faster(locs):
    '''
    locs like [(65564, 65573), (69036, 71585)]
    return a list of positions
    '''
    result_arrays = []
    for s, e in locs:
        result_arrays.append(np.arange(s+1, e+1))  # Create a range for each (s,e) pair
    
    return list(np.concatenate(result_arrays))  # Concatenate all arrays into a single NumPy array


def translateCDSplus(transcript_id, df_transcript2):
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


def getMutations(df_transcript2, df_mutations, transcript_id):
    '''
    given a row in df_transcript2 of transcript_id, return the index of mutations where it locates inside genomicStart and genomicEnd
    Note the genomicStart should be add 1, as it is index starting from 0
    '''
    r = df_transcript2.loc[transcript_id]
    tdf = df_mutations
    genomicLocs = r['genomicLocs']
    results = []
    for genomicStart, genomicEnd in genomicLocs:
        genomicStart = genomicStart + 1
        tdf1 = tdf[(tdf['pos'] >= genomicStart) & (tdf['pos'] <= genomicEnd) & (tdf['pos_end'] >= genomicStart) & (tdf['pos_end'] <= genomicEnd)]
        results += list(tdf1.index)
    return results

def get_df_transcript2(file_gtf, file_protein, file_genome, cpu_counts, datatype, chromosome, df_mutations = None):
    '''
    get df_transcripts2, add AA_seq, add CDSplus and other notes
    for some reason, the strand of sequences in the gtf file is not accurate. need to calculate
    '''
    # get df_gtf
    df_gtf = pd.read_csv(file_gtf, sep='\t',header=None, comment='#')
    df_gtf.columns = ['seqname','source','feature','start','end','score','strand','frame','protein_id']
    # change df_gtf seqname to str
    df_gtf['seqname'] = df_gtf['seqname'].astype(str)
    
    # get df_protein
    df_protein = parse_proteins(file_protein)
    
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
        # 20220505 fix a bug. Need to correct strand in df_gtf
        df_gtf['strand'] = df_gtf['protein_id'].map(strand_n.to_dict())


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
    df_transcript2['seqname'] = chromosome

    # add column "AA_translate"
    pool = Pool(cpu_counts)
    df_transcript2['AA_translate'] = pool.starmap(translateCDSplus, [[i, df_transcript2] for i in df_transcript2.index])
    pool.close()
    pool.join()
    # change 'frame' to int
    df_transcript2['frame'] = df_transcript2['frame'].astype(int)
    
    if df_mutations is None:
        return df_transcript2
    # assign mutations to each transcript. 
    # add column 'mutations' to df_transcript2, which stores index on the mutations in df_mutations
    # add column 'mutations', with list of mutation index in df_mutations
    if df_mutations is not None:
        pool = Pool(cpu_counts)
        results = pool.starmap(getMutations, [[df_transcript2, df_mutations, i] for i in df_transcript2.index])
        pool.close()
        pool.join()
        df_transcript2['mutations'] = results
    
        
    # for RefSeq, sometimes the GTF annotation does not agree with the protein sequence, which means the len(CDS)/len(protein) != 3, find those cases and do not change the proteins
    if datatype == 'RefSeq':
        tdf_special = df_transcript2.loc[~df_transcript2.apply(checkIfAAtranslatedFromGenome,axis=1)]
        if tdf_special.shape[0] > 0:
            print('chromosome', chromosome, list(tdf_special.index), 'not translated from the CDS sequences in the genome. do not change')
            for transcript_id in tdf_special.index:
                df_transcript2.at[transcript_id, 'mutations'] = []
    
    # split df_transcript2 to with mutations and no_mutations
    df_transcript_noMut = df_transcript2[df_transcript2['mutations'].apply(lambda x:len(x)==0)]
    df_transcript3 = df_transcript2[~df_transcript2.index.isin(set(df_transcript_noMut.index))].copy()
    
    print('transcripts with no mutations:', df_transcript_noMut.shape[0])
    print('transcripts with mutations in CDSplus:', df_transcript3.shape[0])

    return df_transcript3


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


def getMut_helper(mutations, strand, df_mutations):
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


def create_df_CDSplus_for_transcript_id(r):
    '''
    r is a row in df_transcript2
    create a dataframe for CDSplus, with chr, pos, ref, alt, variant_id
    '''
    chromosome = r['seqname']
    transcript_id = r.name
    locs = r['genomicLocs']
    locs = getPositionsFromLocs_faster(locs)
    CDSplus = Seq(r['CDSplus'])
    AA_seq = r['AA_seq'].replace('*','_')# in the middle of some sequences, the stop codon is marked as *
    AA_ori = AA_seq
    AA_translate = r['AA_translate']
    frame = int(r['frame'])
    strand = r['strand']
    # mutations = r['mutations']
    # tdf_m = self.getMut(mutations, strand)
    
    # sometimes the AA sequences for example ENSP00000466819, is longer than the CDSplus. skip.
    if len(AA_seq.strip('X'))*3 > len(CDSplus):
        print(transcript_id, 'input protein sequences cannot be translated from the CDS sequence in gtf annotation.')

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
        # tdc_result['nonStandardStopCodon'] = 1
        nonStandardStopCodon = 1
        CDSplus = CDSplus[:3*AA_len]
        locs = locs[:3*AA_len]
    else:
        nonStandardStopCodon = 0
    
    df_CDSplus = pd.DataFrame()
    df_CDSplus['locs'] = locs
    df_CDSplus['bases'] = list(CDSplus)
    df_CDSplus['chr'] = chromosome
    df_CDSplus['strand'] = strand
    return df_CDSplus, AA_seq, AA_ori, AA_translate, nonStandardStopCodon




def translateCDSplusWithMut(r, df_mutations):
    '''
    r is a row of df_transcript3
    transcript_id is index in df_transcript3
    translate CDSplus based on CDSplus and frame, and mutations
    '''
    transcript_id = r.name
    chromosome = r['seqname']
    tdc_result = {}
    frame = int(r['frame'])
    strand = r['strand']
    mutations = r['mutations']
    tdf_m = getMut_helper(mutations, strand, df_mutations)
    
    df_CDSplus, AA_seq, AA_ori, AA_translate, nonStandardStopCodon = create_df_CDSplus_for_transcript_id(r)
    tdc_result['nonStandardStopCodon'] = nonStandardStopCodon
    AA_len = len(AA_seq)
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
    
    # if no special codon, translate directly
    if AA_seq in AA_translate:
        # AA_seq may be shorter than AA_translate, because for annotation, the sequence may end before the stop codon
        CDSnew = ''.join(df_CDSplus['new_nt'])
        new_AA = str(Seq(CDSnew[:3*(len(CDSnew) // 3)]).translate().split('*')[0])
        tdc_result['new_AA'] = new_AA
        # no AA changed, return original sequence
        if new_AA == AA_seq:
            if frame != 0:
                tdc_result['new_AA'] = 'X' + tdc_result['new_AA']
            tdc_result['new_AA'] = tdc_result['new_AA'].replace('_','*')
            return tdc_result
        # only AA change
        if len(new_AA) == len(AA_seq):
            t_alt = pd.DataFrame()
            for n in range(len(AA_seq)):
                if AA_seq[n] != new_AA[n]:
                    t_alt.loc[n, 'AA_ref'] = AA_seq[n]
                    t_alt.loc[n, 'AA_alt'] = new_AA[n]
                    t_alt.loc[n, 'AA_index'] = n + 1
                    t_alt.loc[n, 'variants'] = ','.join(df_CDSplus.iloc[n*3:(n+1)*3]['variant_id'].dropna().drop_duplicates())
            if frame !=0: 
                t_alt['AA_index'] = t_alt['AA_index'] + 1
            tdc_result['n_variant_AA'] = t_alt.shape[0]
            tdc_result['variant_AA'] = ';'.join(t_alt.apply(lambda x:x['AA_ref'] + str(int(x['AA_index'])) + x['AA_alt'] + '({})'.format(x['variants']), axis=1))
            if frame != 0:
                tdc_result['new_AA'] = 'X' + tdc_result['new_AA']
            tdc_result['new_AA'] = tdc_result['new_AA'].replace('_','*')
            return tdc_result
    else:
        if 'U' in AA_seq:
            print(f'{transcript_id} with selenocysteine, special case')
        else:
            print(f'warning! protein {transcript_id}, the original provided protein sequence is {AA_seq} and translated protein sequence from the GTF annotation is {AA_translate}')
            
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
    if frame !=0: 
        codons_ref['AA_index'] = codons_ref['AA_index'] + 1
    # add ref AA to codon_alt
    codons_alt = codons_alt.merge(codons_ref, left_on=['chr','strand','codon1'], right_on=['chr','strand','codon1'], how='left')
    codons_alt.columns = ['chr','strand','codon1', 'codon_alt','variants', 'codon_ref', 'variants_ref','AA_ref','AA_index']
    # 20220510: start codon frame shift complex case, ignore for now, then AA_index column will be all "NaN"
    if all(codons_alt['AA_index'].isnull()):
        print('frame shift, complex case, no change for', transcript_id)
        return {}

    codons_alt.loc[codons_alt.duplicated('AA_index'),['codon_ref','AA_ref','AA_index']] = np.NaN # only use each AA_ref Once
    codons_alt['AA_alt'] = codons_alt.apply(lambda x:x['AA_ref'] if x['codon_alt'] == x['codon_ref'] else str(Seq(x['codon_alt']).translate()), axis=1)
    # codons_alt['AA_index'] = codons_alt['AA_index'].fillna(method='ffill')
    codons_alt['AA_index'] = codons_alt['AA_index'].ffill()
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
    t_alt = t_alt[((t_alt['AA_alt'] == '*') & (t_alt['AA_ref'] != '*')) & ((t_alt['variants'] != '') | t_alt['variants_ref'].notnull())]
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
    # print(time.time() - t0)
    if not frameChange:
        if n_variant_AA + n_deletion_AA + n_insertion_AA == 0 and (not stopGain):
            if AA_ori in new_AA:#new_AA is longer than AA_ori at the end of the sequence, stop loss
                new_AA = new_AA #AA_ori
            elif new_AA in AA_translate:
                print(transcript_id, 'AA_translate longer than provided AA sequence, the gtf file not consistent with the protein file')
            else:
                print(transcript_id,'something wrong, likely to be mutation in start codon?')
                print('input AA', AA_ori)
                print('translated AA', new_AA)
    
    tdc_result['new_AA'] = new_AA.replace('_','*')
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

def translateCDSplusWithMut2(r, df_mutations):
    transcript_id = r.name
    try:
        return translateCDSplusWithMut(r, df_mutations)
    except Exception as e:
        print(f'{transcript_id} cannot be processed properly, please check Error: {e}')
        return {}


def save_mutation_and_proteins(df_transcript3, outprefix):
    '''
    save results based on outprefix
    '''
    # save mutation annotation
    columns_keep = ['protein_id_fasta', 'seqname', 'strand','frameChange','stopGain', 'AA_stopGain', 'stopLoss', 'stopLoss_pos', 'nonStandardStopCodon', 'n_variant_AA', 'n_deletion_AA', 'n_insertion_AA', 'variant_AA', 'insertion_AA', 'deletion_AA', 'len_ref_AA', 'len_alt_AA']
    columns_keep = [e for e in columns_keep if e in df_transcript3.columns]
    df_sum_mutations = df_transcript3[(df_transcript3['AA_seq'] != df_transcript3['new_AA']) & (pd.notnull(df_transcript3['new_AA']))][columns_keep]
    
    outfilename = outprefix +'.aa_mutations.csv'
    if not os.path.exists(os.path.dirname(outfilename)):
        os.makedirs(os.path.dirname(outfilename))
    df_sum_mutations.to_csv(outfilename, sep='\t')
    print('number of proteins with AA change:', df_sum_mutations.shape[0])
    
    # save proteins
    fout = open(outprefix + '.mutated_protein.fa','w')
    for protein_id, r in df_transcript3.iterrows():
        if pd.notnull(r['new_AA']) and r['new_AA'] != r['AA_seq']:
            fout.write('>{}\tchanged\n{}\n'.format(r['protein_description'], r['new_AA']))
    fout.close()

class PerChrom(object):
    """
    Define an instance of PerChromosome analysis
    do perGeno analysis for data from the same chromosome
    """
    def __init__(
                    self,
                    file_genome,
                    file_gtf,
                    file_mutations,
                    file_protein,
                    threads,
                    outprefix,
                    datatype,
                    chromosome
                ):
        # variable from input
        self.file_genome = file_genome # genome file location
        self.file_gtf = file_gtf # gtf file location
        self.file_mutations = file_mutations # mutation file location
        self.file_protein = file_protein # protein file location
        self.threads = threads # threads to use when running program
        self.outprefix = outprefix # where to store the results
        self.datatype = datatype # input datatype, could be GENCODE_GTF, GENCODE_GFF3, RefSeq or gtf
        self.chromosome = chromosome # chromosome name
        
        # dataframe to store the mutation information
        if self.file_mutations:
            self.df_mutations = parse_mutation(file_mutations=self.file_mutations, chromosome=self.chromosome)
        else:
            print(self.chromosome, 'No mutation file provided, will not do mutation analysis')
            self.df_mutations = None
        # dataframe to store all transcript info
        self.df_transcript2 = get_df_transcript2(file_gtf=self.file_gtf, file_protein=self.file_protein, file_genome=self.file_genome,cpu_counts=self.threads,datatype=self.datatype, chromosome=self.chromosome, df_mutations=self.df_mutations)


    def run_perChrom(self, save_results = True):
        '''run perChrom
        '''
        cpu_counts = self.threads
        df_transcript3 = self.df_transcript2
        datatype = self.datatype
        outprefix = self.outprefix
        chromosome = self.chromosome
        df_mutations = self.df_mutations
        
        if df_transcript3.shape[0] == 0:
            print('No protein sequences to change for chromosome', chromosome)
            return df_transcript3

        pool = Pool(cpu_counts)
        results = pool.starmap(translateCDSplusWithMut2, [[r, df_mutations] for _,r in df_transcript3.iterrows()])
        pool.close()
        pool.join()
        tdf = pd.DataFrame(results)
        for col in tdf.columns:
            df_transcript3[col] = list(tdf[col])
        
        # get sequence length of ref_AA and alt_AA
        df_transcript3['len_ref_AA'] = df_transcript3['AA_seq'].str.len()
        if 'new_AA' in df_transcript3.columns:
            df_transcript3['len_alt_AA'] = df_transcript3['new_AA'].str.len()
        
        if save_results:
            save_mutation_and_proteins(df_transcript3, outprefix)
        else:
            return df_transcript3



description = '''output a new reference protein set by with the variants data for each chromosome. The input files were generated by PrecisionProDB_core.
'''

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-g','--genome', help = 'the reference genome sequence in fasta format. It can be a gzip file', required=True)
    parser.add_argument('-p','--protein', help = 'protein sequences in fasta format. It can be a gzip file.', default=None)
    parser.add_argument('-f', '--gtf', help='gtf file with CDS and exon annotations', required=True)
    parser.add_argument('-m', '--mutations', help='a file stores the variants', required=True)
    parser.add_argument('-t', '--threads', help='number of threads/CPUs to run the program. default, use all CPUs available', type=int, default=os.cpu_count())
    parser.add_argument('-o', '--out', help='output prefix. two file will be output. One is the annotation for mutated transcripts, one is the protein sequences. {out}.aa_mutations.csv, {out}.mutated_protein.fa. {out}.mutated_protein.fa only includes the proteins AA changes. ''')
    parser.add_argument('-c', '--chromosome', help='''chromosome name/id, default="chr1" ''', default='chr1', type=str)
    parser.add_argument('-a', '--datatype', help='''input datatype, could be GENCODE_GTF, GENCODE_GFF3, RefSeq, Ensembl_GTF or gtf. default "gtf". Ensembl_GFF3 is not supported. ''', default='gtf', type=str, choices=['GENCODE_GTF', 'GENCODE_GFF3','RefSeq','Ensembl_GTF','gtf'])
    
    f = parser.parse_args()
    
    file_genome = f.genome
    file_protein = f.protein
    file_gtf = f.gtf
    file_mutations = f.mutations
    outprefix = f.out
    threads = f.threads
    chromosome = f.chromosome
    datatype = f.datatype
    
    perchrom = PerChrom(file_genome = file_genome,
                    file_gtf = file_gtf,
                    file_mutations = file_mutations,
                    file_protein = file_protein,
                    threads = threads,
                    outprefix = outprefix,
                    datatype = datatype,
                    chromosome = chromosome)
    print('run perChrom for chromosome', chromosome)
    perchrom.run_perChrom()

