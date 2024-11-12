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

# Global variables
con = None
df_mutations = None

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



def process_mutation_row(row_index):
    global con, df_mutations
    row = df_mutations.loc[row_index]
    chromosome, pos = row['chr'], row['pos']
    protein_ids = buildSqlite.get_protein_id_from_genomicLocs(con, chromosome, pos)
    return [row.name, protein_ids]

def get_protein_id_from_df_mutations(cpu_counts=10):
    '''
    get protein_id from df_mutations
    df_mutation is a dataframe with the mutations
    return a dictionary, with protein_id as keys, and list of df_mutation index as values
    '''
    global df_mutations, con
    pool = Pool(cpu_counts)
    results = pool.map(process_mutation_row, list(df_mutations.index))
    pool.close()
    pool.join()
    
    combined_results = {}
    for indices, protein_ids in results:
        for protein_id in protein_ids:
            if protein_id not in combined_results:
                combined_results[protein_id] = []
            combined_results[protein_id].append(indices)
    
    return combined_results


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
                    chromosome
                ):
        self.file_sqlite = file_sqlite # genome file location
        self.file_mutations = file_mutations # mutation file location
        self.threads = threads # threads to use when running program
        self.outprefix = outprefix # where to store the results
        self.datatype = datatype # input datatype, could be GENCODE_GTF, GENCODE_GFF3, RefSeq or gtf
        self.chromosome = chromosome # chromosome name
        
        # dataframe to store the mutation information
        global df_mutations
        
        if isinstance(self.file_mutations, pd.DataFrame):
            df_mutations = self.file_mutations
            self.df_mutations = df_mutations
            # print(df_mutations)
        elif self.file_mutations:
            df_mutations = perChrom.parse_mutation(file_mutations=self.file_mutations, chromosome=self.chromosome)
            # print(df_mutations.head())
            self.df_mutations = df_mutations
        else:
            print(self.chromosome, 'No mutation file provided, will not do mutation analysis')
        
        global con
        if con is None:
            con = buildSqlite.get_connection(self.file_sqlite)


    def run_perChrom(self, save_results = True):
        '''run perChrom
        '''
        cpu_counts = self.threads
        datatype = self.datatype
        outprefix = self.outprefix
        chromosome = self.chromosome
        df_mutations = self.df_mutations
        global con

        # get transcripts that will be used based on the self_mutations
        dc_transcript2mutations = get_protein_id_from_df_mutations(cpu_counts=cpu_counts)

        # get df_transcript2 based on dc_transcript2mutations
        protein_ids = tuple(dc_transcript2mutations.keys())
        query = f"SELECT * FROM protein_description WHERE protein_id IN {protein_ids}"
        df_transcript2 = pd.read_sql_query(query, con)        # assign mutations to each transcript. 
        
        # update df_transcript2
        df_transcript2 = df_transcript2.set_index('protein_id')
        df_transcript2['genomicLocs'] = df_transcript2['genomicLocs'].apply(eval)
        # add column 'mutations', with list of mutation index in df_mutations
        pool = Pool(cpu_counts)
        results = pool.starmap(perChrom.getMutations, [[df_transcript2, df_mutations, i] for i in df_transcript2.index])
        pool.close()
        pool.join()
        df_transcript2['mutations'] = results
        
        # for RefSeq, sometimes the GTF annotation does not agree with the protein sequence, which means the len(CDS)/len(protein) != 3, find those cases and do not change the proteins
        if datatype == 'RefSeq':
            tdf_special = df_transcript2.loc[~df_transcript2.apply(perChrom.checkIfAAtranslatedFromGenome,axis=1)]
            if tdf_special.shape[0] > 0:
                print('chromosome', chromosome, list(tdf_special.index), 'not translated from the CDS sequences in the genome. do not change')
                for transcript_id in tdf_special.index:
                    df_transcript2.at[transcript_id, 'mutations'] = []

        df_transcript3 = df_transcript2
    
        if df_transcript3.shape[0] == 0:
            print('No protein sequences to change for chromosome', chromosome)
            return df_transcript3


        pool = Pool(cpu_counts)
        results = pool.starmap(perChrom.translateCDSplusWithMut2, [[r, df_mutations] for _,r in df_transcript3.iterrows()])
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
            perChrom.save_mutation_and_proteins(df_transcript3, outprefix)
        else:
            return df_transcript3



description = '''output a new reference protein set by with the variants data for each chromosome. The input files were generated by PrecisionProDB_core.
'''

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-q', '--sqlite', help='sqlite file file with CDS and exon annotations', required=True)
    parser.add_argument('-m', '--mutations', help='a file stores the variants, or a string looks like "1-55051215-G-GA" or "chr1-55051215-G-GA" ', required=True)
    parser.add_argument('-t', '--threads', help='number of threads/CPUs to run the program. default, use all CPUs available', type=int, default=os.cpu_count())
    parser.add_argument('-o', '--out', help='output prefix. two file will be output. One is the annotation for mutated transcripts, one is the protein sequences. {out}.aa_mutations.csv, {out}.mutated_protein.fa. default "perChrom" ', default = "perChrom")
    parser.add_argument('-c', '--chromosome', help='''chromosome name/id, default="chr1" ''', default='chr1', type=str)
    parser.add_argument('-a', '--datatype', help='''input datatype, could be GENCODE_GTF, GENCODE_GFF3, RefSeq, Ensembl_GTF or gtf. default "gtf". Ensembl_GFF3 is not supported. ''', default='gtf', type=str, choices=['GENCODE_GTF', 'GENCODE_GFF3','RefSeq','Ensembl_GTF','gtf'])
    
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
    
    con = buildSqlite.get_connection(file_sqlite)
    df_mutations = perChrom.parse_mutation(file_mutations=file_mutations, chromosome=chromosome)

    perchrom_sqlite = PerChrom_sqlite(file_sqlite = file_sqlite,
                    file_mutations = file_mutations,
                    threads = threads,
                    outprefix = outprefix,
                    datatype = datatype,
                    chromosome = chromosome)
    print('run perChrom for chromosome', chromosome)
    perchrom_sqlite.run_perChrom()

