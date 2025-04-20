import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import gzip
import numpy as np
from multiprocessing import Pool
import os
import time
import perChrom
import re
import buildSqlite
import perChrom
import sqlite3
import perChrom
try:
    import tqdm
except:
    print('Cannot import tqdm. Will not use tqdm.')
    tqdm = False
# # Global variables
# con = None
# df_mutations = None

# code below for testing the the program
TEST = False
# TEST = True
if TEST:
    file_mutations = '/XCLabServer002_fastIO/examples/gnomAD.variant.txt.gz'
    file_mutations = 'chr1-973858-G-C'
    threads = 10
    outprefix = '/XCLabServer002_fastIO/examples/GENCODE/GENCODE.tsv_temp/chr1'
    datatype = 'GENCODE_GTF'
    file_sqlite = '/XCLabServer002_fastIO/examples/GENCODE/GENCODE.tsv.sqlite'
    chromosome = 'chr1'
    con = buildSqlite.get_connection(file_sqlite)
    df_mutations = perChrom.parse_mutation(file_mutations=file_mutations, chromosome=chromosome)



# def process_mutation_row(row_index):
#     # global con, df_mutations
#     row = df_mutations.loc[row_index]
#     chromosome, pos = row['chr'], row['pos']
#     protein_ids = buildSqlite.get_protein_id_from_genomicLocs(con, chromosome, pos)
#     return [row.name, protein_ids]

def is_valid(value):
    """
    Check if a value is not None, '', 0, '0', pd.NaT, pd.NaN, or np.nan.
    
    Parameters:
    value: The value to check.
    
    Returns:
    bool: True if the value is valid, False otherwise.
    """
    if pd.isna(value) or value in [None, '', 0, '0', [], 'False', False, '.', 'NaN', '*', 'N/A', 'NA', 'nan', 'null', 'Null', 'NULL', 'None', 'none','-', ' ', '[]']:
        return False

    if value in [1, '1', True, 'True', 2, '2', 3, '3', 4, '4', 5, '5', 6, '6', 7, '7', 8, '8', 9, '9']:
        return True

    return True


def get_protein_id_from_df_mutations(df_mutations, file_sqlite, cpu_counts=10):
    '''
    get protein_id from df_mutations
    df_mutation is a dataframe with the mutations
    return a dictionary, with protein_id as keys, and list of df_mutation index as values
    '''
    con = buildSqlite.get_connection(file_sqlite)
    if df_mutations.shape[0] < 1000:
        threads = None
    else:
        threads = cpu_counts

    protein_ids = buildSqlite.get_protein_id_from_genomicLocs(con, list(df_mutations['chr']), list(df_mutations['pos']), list(df_mutations['pos_end']), threads = threads)
    
    results = [[i,j] for i,j in zip(df_mutations.index, protein_ids)]
    
    combined_results = {}
    for indices, protein_ids in results:
        for protein_id in protein_ids:
            if protein_id not in combined_results:
                combined_results[protein_id] = []
            combined_results[protein_id].append(indices)
    
    return combined_results

def convert_df_transcript2_to_df_transcript3_helper(protein_id, df_transcript2, df_mutations, individual, kargs):
    '''
    protein_id is a protein_id in df_transcript2
    for each protein_id, get the mutations in each individual.
    return a list of tuple, each tuple is (tuple of variant_index, individuals with this variant_index joined by ',')
    
    If shape and file_memmap are provided in kargs, it means we're dealing with an extra large mutation file
    where individual data is stored in a memory-mapped file instead of in the DataFrame.
    '''
    mutations = df_transcript2.loc[protein_id]['mutations']
    
    # Check if we're using memory-mapped file for individual data
    if 'shape' in kargs and 'file_memmap' in kargs:
        # Using memory-mapped file for individual data
        shape = kargs['shape']
        file_memmap = kargs['file_memmap']
        
        # Open the memory-mapped file in read mode
        mmap = np.memmap(file_memmap, dtype='int8', mode='r', shape=shape)
        
        # Get the indices of mutations for this protein
        mutation_indices = mutations
        
        tdc = {}
        for i, sample in enumerate(individual):
            # For each sample, check which mutations are valid (value is 1)
            valid_mutations = []
            for idx in mutation_indices:
                if idx < mmap.shape[0] and i < mmap.shape[1] and mmap[idx, i] == 1:
                    valid_mutations.append(idx)
            
            variant_index = tuple(valid_mutations)
            if len(variant_index) > 0:
                if variant_index not in tdc:
                    tdc[variant_index] = []
                tdc[variant_index].append(sample)
    else:
        # Using regular DataFrame for individual data
        tdf_m = df_mutations.loc[mutations]
        tdc = {}
        for sample in individual:
            if sample not in tdf_m.columns:
                print('warning:', sample, 'not in mutation columns and is ignored! unexpected error may happen!')
                continue
            tdf_m_single = tdf_m[tdf_m[sample].apply(is_valid)]
            variant_index = tuple(tdf_m_single.index)
            if len(variant_index) > 0:
                if variant_index not in tdc:
                    tdc[variant_index] = []
                tdc[variant_index].append(sample)
    
    tdc = {k:','.join(v) for k,v in tdc.items()}
    return list(tdc.items())

def convert_df_transcript2_to_df_transcript3(df_transcript2, df_mutations, individual = None, cpu_counts=1, **kargs):
    '''
    convert df_transcript2 to df_transcript3
    if individual is None or individual is '', return df_transcript2
    else, group mutations in df_transcript2 by individual. rename index by adding __X, and add column 'individual'
    individual is a list of samples in df_mutations
    
    If shape and file_memmap are provided in kargs, it means we're dealing with an extra large mutation file
    where individual data is stored in a memory-mapped file instead of in the DataFrame.
    '''
    if individual is None or individual == '' or individual == 'None' or individual == [] or individual == 'ALL_VARIANTS':
        return df_transcript2
    
    if isinstance(individual, list):
        individual = individual
    elif individual == 'ALL_SAMPLES':
        individual = [i for i in df_mutations.columns if i not in ['chr', 'pos', '', 'ref', 'alt', 'pos_end']]
        if len(individual) < 100:
            print('warning: individual is ALL_SAMPLES, all columns in file_mutations other than chr pos ref alt were used as individuals, which may cause problem. individuals were set as:', individual)
        else:
            print('warning: individual is ALL_SAMPLES, all columns in file_mutations other than chr pos ref alt were used as individuals, which may cause problem. {} individuals were used'.format(len(individual)))
    elif ',' in individual:
        individual = individual.split(',')
    else:
        individual = [individual]
    
    if cpu_counts == 1:
        results = [convert_df_transcript2_to_df_transcript3_helper(protein_id, df_transcript2, df_mutations, individual, kargs) for protein_id in df_transcript2.index]
    else:
        pool = Pool(cpu_counts)
        results = pool.starmap(convert_df_transcript2_to_df_transcript3_helper, [(protein_id, df_transcript2, df_mutations, individual, kargs) for protein_id in df_transcript2.index])
        pool.close()
        pool.join()
    
    df_transcript3 = df_transcript2.copy()
    df_transcript3['mutations_anno'] = results
    df_transcript3 = df_transcript3.explode('mutations_anno')
    df_transcript3['mutations'] = df_transcript3['mutations_anno'].str[0]
    df_transcript3['individual'] = df_transcript3['mutations_anno'].str[1]
    df_transcript3 = df_transcript3.drop('mutations_anno', axis=1)
    df_transcript3 = df_transcript3[df_transcript3['individual'].apply(is_valid)]
    df_transcript3['mutations'] = df_transcript3['mutations'].apply(list)
    return df_transcript3

def save_mutation_and_proteins(df_transcript3, outprefix):
    '''
    save results based on outprefix
    '''
    # save mutation annotation
    columns_keep = ['protein_id_fasta', 'seqname', 'strand','frameChange','stopGain', 'AA_stopGain', 'stopLoss', 'stopLoss_pos', 'nonStandardStopCodon', 'n_variant_AA', 'n_deletion_AA', 'n_insertion_AA', 'variant_AA', 'insertion_AA', 'deletion_AA', 'len_ref_AA', 'len_alt_AA','individual','new_AA','AA_seq']
    columns_keep = [e for e in columns_keep if e in df_transcript3.columns]
    if df_transcript3.shape[0] == 0:
        print('no protein with AA change')
        return None
    df_sum_mutations = df_transcript3[(df_transcript3['AA_seq'] != df_transcript3['new_AA']) & (pd.notnull(df_transcript3['new_AA']))][columns_keep]
    
    df_sum_mutations = df_sum_mutations.reset_index()
    df_sum_mutations = df_sum_mutations.groupby([i for i in df_sum_mutations.columns if i != 'individual'] ,dropna=False)['individual'].apply(lambda x:','.join(x)).reset_index()
    df_sum_mutations['protein_id_fasta_nth'] = df_sum_mutations.groupby('protein_id_fasta').cumcount()+1
    df_sum_mutations['protein_id_fasta'] = df_sum_mutations.apply(lambda x: '{}__{}'.format(x['protein_id_fasta'], x['protein_id_fasta_nth']), axis=1)
    outfilename = outprefix +'.aa_mutations.csv'
    if not os.path.exists(os.path.dirname(outfilename)):
        os.makedirs(os.path.dirname(outfilename))
    df_sum_mutations[[i for i in df_sum_mutations.columns if i not in ['new_AA','AA_seq', 'protein_id_fasta_nth']]].to_csv(outfilename, sep='\t',index=None)
    print('{} proteins with AA change and generated {} mutated proteins'.format(df_sum_mutations['protein_id'].nunique(), df_sum_mutations['protein_id_fasta'].nunique()))
    
    
    # save proteins
    fout = open(outprefix + '.mutated_protein.fa','w')
    for protein_id, r in df_sum_mutations.iterrows():
        if pd.notnull(r['new_AA']) and r['new_AA'] != r['AA_seq']:
            fout.write('>{}\tchanged\n{}\n'.format(r['protein_id_fasta'], r['new_AA']))
    fout.close()
# Define the wrapper function needed for imap
def translate_wrapper(args):
    """Unpacks arguments for translateCDSplusWithMut2"""
    row_series, mutations_df = args
    # Assuming perChrom.translateCDSplusWithMut2 exists and works like this
    return perChrom.translateCDSplusWithMut2(row_series, mutations_df)


class PerChrom_sqlite(object):
    """
    Define an instance of PerChromosome analysis
    do perGeno analysis for data from the same chromosome
    """
    def __init__(
                    self,
                    file_sqlite,
                    file_mutations,
                    threads,
                    outprefix,
                    datatype,
                    chromosome,
                    individual = None
                ):
        self.file_sqlite = file_sqlite # genome file location
        self.file_mutations = file_mutations # mutation file location
        self.threads = threads # threads to use when running program
        self.outprefix = outprefix # where to store the results
        self.datatype = datatype # input datatype, could be GENCODE_GTF, GENCODE_GFF3, RefSeq or gtf
        self.chromosome = chromosome # chromosome name
        self.extra_large_file_mutation = False # whether the mutation file is larger than 1G
        
        # if file_mutation is larger than 1G, only read ['chr', 'pos', 'ref', 'alt']
        if isinstance(self.file_mutations, str):
            if os.path.exists(self.file_mutations):
                if os.path.getsize(self.file_mutations) > 100000000:
                    self.extra_large_file_mutation = True
        
        if self.extra_large_file_mutation:
            df_mutations = perChrom.readExtraLargeMutationFile(file_mutations, nrows=2)
        elif isinstance(self.file_mutations, pd.DataFrame):
            df_mutations = self.file_mutations
            self.df_mutations = df_mutations
            # print(df_mutations)
        elif self.file_mutations:
            df_mutations = perChrom.parse_mutation(file_mutations=self.file_mutations, chromosome=self.chromosome)
            # print(df_mutations.head())
            self.df_mutations = df_mutations
        else:
            print(self.chromosome, 'No mutation file provided, will not do mutation analysis')
            self.df_mutations = pd.DataFrame(columns=['chr', 'pos', 'ref', 'alt'])
        
        self.con = buildSqlite.get_connection(self.file_sqlite)
        
        if individual is None or individual == '' or individual == 'None' or individual == [] or individual == 'ALL_VARIANTS':
            individual = None
        elif isinstance(individual, list):
            individual = individual
        elif individual == 'ALL_SAMPLES':
            individual = [i for i in df_mutations.columns if i not in ['chr', 'pos', '', 'ref', 'alt', 'pos_end']]
            if len(individual) < 100:
                print('warning: individual is ALL_SAMPLES, all columns in file_mutations other than chr pos ref alt were used as individuals, which may cause problem. individuals were set as:', individual)
            else:
                print('warning: individual is ALL_SAMPLES, all columns in file_mutations other than chr pos ref alt were used as individuals, which may cause problem. {} individuals were used'.format(len(individual)))

        elif ',' in individual:
            individual = individual.split(',')
        elif individual in df_mutations.columns:
            individual = [individual]
        else:
            individual = None
            print('warning: individual is not in file_mutations, will not be used for do mutation analysis')
        
        self.individual = individual
        
        if self.extra_large_file_mutation:
            self.df_mutations = perChrom.parse_mutation(file_mutations, columns_to_include=['chr', 'pos', 'ref', 'alt'])
            from vcf2mutation import tsv2memmap
            shape = (self.df_mutations.shape[0], len(self.individual))
            tsv2memmap(file_mutations, individuals = self.individual, memmap_file=file_mutations + '.memmap')
            self.shape = shape
            self.file_memmap = file_mutations + '.memmap'
        
        if chromosome:
            chromosome_with_genomicLocs = buildSqlite.get_genomicLocs_chromosomes(file_sqlite)
            if chromosome in chromosome_with_genomicLocs:
                self.df_mutations['chr'] = chromosome
            elif chromosome.startswith('chr') and chromosome[3:] in chromosome_with_genomicLocs:
                self.df_mutations['chr'] = chromosome[3:]
            elif 'chr' + chromosome in chromosome_with_genomicLocs:
                self.df_mutations['chr'] = 'chr' + chromosome
            else:
                print(f'warning! chromosome {chromosome} not in file_sqlite')



    def run_perChrom(self, save_results = True):
        '''run perChrom
        '''
        cpu_counts = self.threads
        datatype = self.datatype
        outprefix = self.outprefix
        chromosome = self.chromosome
        df_mutations = self.df_mutations
        individual = self.individual
        # global con

        # get transcripts that will be used based on the self_mutations
        dc_transcript2mutations = get_protein_id_from_df_mutations(df_mutations, file_sqlite=self.file_sqlite, cpu_counts=cpu_counts)

        # get df_transcript2 based on dc_transcript2mutations
        protein_ids = tuple(dc_transcript2mutations.keys())
        if len(protein_ids) == 1:
            protein_ids = '("{}")'.format(protein_ids[0])
        query = f"SELECT * FROM protein_description WHERE protein_id IN {protein_ids}"
        df_transcript2 = pd.read_sql_query(query, self.con)        # assign mutations to each transcript. 
        if df_transcript2.shape[0] == 0:
            print('No protein sequences to change for chromosome', chromosome)
            return df_transcript2
        # update df_transcript2
        df_transcript2 = df_transcript2.set_index('protein_id')
        df_transcript2['genomicLocs'] = df_transcript2['genomicLocs'].apply(eval)
        # add column 'mutations', with list of mutation index in df_mutations
        # pool = Pool(cpu_counts)
        # results = pool.starmap(perChrom.getMutations, [[df_transcript2, df_mutations, i] for i in df_transcript2.index])
        # pool.close()
        # pool.join()
        df_transcript2['mutations'] = [dc_transcript2mutations[i] for i in df_transcript2.index]
        df_transcript2 = df_transcript2[df_transcript2['mutations'].str.len() > 0]
        
        # for RefSeq, sometimes the GTF annotation does not agree with the protein sequence, which means the len(CDS)/len(protein) != 3, find those cases and do not change the proteins
        if datatype == 'RefSeq':
            tdf_special = df_transcript2.loc[~df_transcript2.apply(perChrom.checkIfAAtranslatedFromGenome,axis=1)]
            if tdf_special.shape[0] > 0:
                print('chromosome', chromosome, list(tdf_special.index), 'not translated from the CDS sequences in the genome. do not change')
                for transcript_id in tdf_special.index:
                    df_transcript2.at[transcript_id, 'mutations'] = []
                df_transcript2 = df_transcript2[~df_transcript2.index.isin(set(tdf_special.index))]

        if not self.extra_large_file_mutation:
            df_transcript3 = convert_df_transcript2_to_df_transcript3(df_transcript2, df_mutations, individual = individual, cpu_counts=cpu_counts)
        else:
            df_transcript3 = convert_df_transcript2_to_df_transcript3(df_transcript2, df_mutations, individual = individual, cpu_counts=cpu_counts, shape=self.shape, file_memmap = self.file_memmap)
    
        if df_transcript3.shape[0] == 0:
            print('No protein sequences to change for chromosome', chromosome)
            return df_transcript3

        if cpu_counts > 1:
            pool = Pool(cpu_counts)
            starmap_args = [[r, df_mutations] for _,r in df_transcript3.iterrows()]
            
            total_tasks = df_transcript3.shape[0]
            chunk_size = min(200, total_tasks // cpu_counts // 4)
            
            # results = pool.starmap(perChrom.translateCDSplusWithMut2, starmap_args, chunksize=100)
            imap_results = pool.imap(translate_wrapper, starmap_args, chunksize=chunk_size)
            if tqdm:
                results = list(tqdm.tqdm(imap_results, total=total_tasks, desc=f"Translating {chromosome}"))
            else:
                results = list(imap_results)
            pool.close()
            pool.join()
        else:
            results = [perChrom.translateCDSplusWithMut2(r, df_mutations) for _,r in df_transcript3.iterrows()]
        tdf = pd.DataFrame(results)
        for col in tdf.columns:
            df_transcript3[col] = list(tdf[col])
        
        # get sequence length of ref_AA and alt_AA
        df_transcript3['len_ref_AA'] = df_transcript3['AA_seq'].str.len()
        if 'new_AA' in df_transcript3.columns:
            df_transcript3['len_alt_AA'] = df_transcript3['new_AA'].str.len()
        
        if save_results:
            if individual is None or individual == '' or individual == 'None' or individual == [] or individual == 'ALL_VARIANTS':
                perChrom.save_mutation_and_proteins(df_transcript3, outprefix)
            else:
                save_mutation_and_proteins(df_transcript3, outprefix)
        else:
            return df_transcript3



description = '''output a new reference protein set by with the variants data for each chromosome. The input files were generated by PrecisionProDB_core.
'''

def main():
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-q', '--sqlite', help='sqlite file file with CDS and exon annotations', required=True)
    parser.add_argument('-m', '--mutations', help='a file stores the variants, or a string looks like "1-55051215-G-GA" or "chr1-55051215-G-GA" ', required=True)
    parser.add_argument('-t', '--threads', help='number of threads/CPUs to run the program. default, use all CPUs available', type=int, default=min(20, os.cpu_count()))
    parser.add_argument('-o', '--out', help='output prefix. two file will be output. One is the annotation for mutated transcripts, one is the protein sequences. {out}.aa_mutations.csv, {out}.mutated_protein.fa. default "perChrom" ', default = "perChrom")
    parser.add_argument('-c', '--chromosome', help='''chromosome name/id, default="chr1" ''', default='chr1', type=str)
    parser.add_argument('-a', '--datatype', help='''input datatype, could be GENCODE_GTF, GENCODE_GFF3, RefSeq, Ensembl_GTF or gtf. default "gtf". Ensembl_GFF3 is not supported. ''', default='gtf', type=str, choices=['GENCODE_GTF', 'GENCODE_GFF3','RefSeq','Ensembl_GTF','gtf'])
    parser.add_argument('-s', '--sample', help='''
                        sample name in the vcf to extract the variant information. default: None, use all variants and do not consider samples.
                        For multiple samples, use "," to join the sample names. For example, "--sample sample1,sample2,sample3".
                        To use all samples, use "--sample ALL_SAMPLES". 
                        To use all variants regardless where the variants from, use "--sample ALL_VARIANTS".
                        ''', default=None)
    
    if TEST:
        f = parser.parse_args(["-q", file_sqlite, "-m", file_mutations, "-t", str(threads), "-o", outprefix, "-c", chromosome, "-a", datatype])
    else:
        f = parser.parse_args()

    
    file_sqlite = f.sqlite
    file_mutations = f.mutations
    outprefix = f.out
    threads = f.threads
    chromosome = f.chromosome
    datatype = f.datatype
    individual = f.sample
    
    # con = buildSqlite.get_connection(file_sqlite)
    # df_mutations = perChrom.parse_mutation(file_mutations=file_mutations, chromosome=chromosome)

    perchrom_sqlite = PerChrom_sqlite(file_sqlite = file_sqlite,
                    file_mutations = file_mutations,
                    threads = threads,
                    outprefix = outprefix,
                    datatype = datatype,
                    chromosome = chromosome,
                    individual = individual)
    print('run perChrom for chromosome', chromosome)
    perchrom_sqlite.run_perChrom()

if __name__ == '__main__':
    main()