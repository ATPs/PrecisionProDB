import pandas as pd
from Bio import SeqIO
import gzip
import os
import time
import pickle
import perChrom
import shutil
import re
import sys
import perChromSqlite
import buildSqlite
from PrecisionProDB_core import PerGeno, get_k_new
import re
import sqlite3

# code below for testing the the program
TEST = False
# TEST = True
if TEST:
    file_genome = '/XCLabServer002_fastIO/examples/RefSeq/RefSeq.genome.fa.gz'
    file_protein = '/XCLabServer002_fastIO/examples/RefSeq/RefSeq.protein.fa.gz'
    file_gtf = '/XCLabServer002_fastIO/examples/RefSeq/RefSeq.gtf.gz'
    datatype = 'RefSeq'
    file_genome = '/XCLabServer002_fastIO/examples/GENCODE/GENCODE.genome.fa.gz'
    file_protein = '/XCLabServer002_fastIO/examples/GENCODE/GENCODE.protein.fa.gz'
    file_gtf = '/XCLabServer002_fastIO/examples/GENCODE/GENCODE.gtf.gz'
    datatype = 'GENCODE_GTF'
    file_mutations = '/XCLabServer002_fastIO/examples/gnomAD.variant.txt.gz'
    file_mutations = '/XCLabServer002_fastIO/examples/celline.vcf.gz'
    threads = 10
    outprefix = '/XCLabServer002_fastIO/examples/GENCODE/GENCODE.tsv'
    protein_keyword = 'auto'
    filter_PASS = True
    individual = None
    chromosome_only = True
    keep_all = True
    file_sqlite = '/XCLabServer002_fastIO/examples/GENCODE/GENCODE.tsv.sqlite'


def runSinglePerChromSqlite(file_sqlite, file_mutations, tempfolder, threads, chromosome, datatype, individual):
    '''
    run PerChrom_sqlite for a single chromosome
    '''
    outprefix = os.path.join(tempfolder, chromosome)
    # print(file_mutations)
    perchrom_sqlite = perChromSqlite.PerChrom_sqlite(file_sqlite = file_sqlite,
                    file_mutations = file_mutations,
                    threads = threads,
                    outprefix = outprefix,
                    datatype = datatype,
                    chromosome = chromosome,
                    individual = individual)
    print('run perchrom_sqlite for chromosome', chromosome)
    # perchrom_sqlite.run_perChrom()
    # print('finished running perchrom_sqlite for chromosome:', chromosome)
    # return chromosome
    try:
        perchrom_sqlite.run_perChrom()
        print('finished running perchrom_sqlite for chromosome:', chromosome)
        return chromosome
    except Exception as e:
        print('cannot run perchrom_sqlite for chromosome', chromosome, 'Proteins will be unchanged. Error message:', e)
        return None


def runPerChomSqlite(file_sqlite, file_mutations, threads, outprefix, protein_keyword, datatype, keep_all, individual, chromosomes_genome, chromosomes_genome_description, file_gtf):
    '''
    '''
    # split file_mutations by chromosome
    pergeno = PerGeno(
        file_genome = '', 
        file_gtf = file_gtf, 
        file_mutations = file_mutations, 
        file_protein='', 
        threads=threads, 
        outprefix=outprefix, 
        datatype=datatype, 
        protein_keyword=protein_keyword, 
        keep_all = keep_all
        )
    
    if individual == 'ALL_VARIANTS':
        individual = ''
    tempfolder = pergeno.tempfolder
    pergeno.chromosomes_genome = chromosomes_genome
    chromosomes_mutation = pergeno.splitMutationByChromosomeLarge(chromosomes_genome_description=chromosomes_genome_description, chromosomes_genome=chromosomes_genome)

    # run runSinglePerChromSqlite
    chromosomes_mutated = [runSinglePerChromSqlite(file_sqlite, f'{tempfolder}/{chromosome}.mutation.tsv', tempfolder, threads, chromosome, datatype, individual) for chromosome in chromosomes_mutation]
    # successful chromosomes
    chromosomes_mutated = [e for e in chromosomes_mutated if e is not None]
    # collect mutation annotations
    files_mutAnno = ['{}/{}.aa_mutations.csv'.format(tempfolder, chromosome) for chromosome in chromosomes_mutated]
    file_mutAnno = outprefix + '.pergeno.aa_mutations.csv'
    df_mutAnno = pd.concat([pd.read_csv(f, sep='\t') for f in files_mutAnno if os.path.exists(f)], ignore_index=True)
    print('total number of proteins with AA mutation:', df_mutAnno.shape[0])
    df_mutAnno.to_csv(file_mutAnno, sep='\t', index=None)

    # collect protein sequences
    files_proteins_changed = ['{}/{}.mutated_protein.fa'.format(tempfolder, chromosome) for chromosome in chromosomes_mutated]
    file_proteins_changed = outprefix + '.pergeno.protein_changed.fa'
    file_proteins_all = outprefix + '.pergeno.protein_all.fa'
    fout_proteins_changed = open(file_proteins_changed, 'w')
    fout_proteins_all = open(file_proteins_all, 'w')
    proteins_changed_ids = []
    for f in files_proteins_changed:
        if not os.path.exists(f):
            continue
        for s in SeqIO.parse(f,'fasta'):
            fout_proteins_all.write('>{}\n{}\n'.format(s.description, str(s.seq)))
            if s.description.endswith('\tchanged'):
                proteins_changed_ids.append(s.id)
                fout_proteins_changed.write('>{}\n{}\n'.format(s.description, str(s.seq)))
        
    proteins_changed_ids = set([i.split('__')[0] for i in proteins_changed_ids])
    conn = buildSqlite.get_connection(file_sqlite)
    df_protein_description = pd.read_sql("SELECT * FROM protein_description", conn)
    df_protein_description = df_protein_description[~ df_protein_description['protein_id_fasta'].isin(proteins_changed_ids)]
    for _, s in df_protein_description.iterrows():
        protein_id, protein_description, protein_id_fasta, AA_seq = s['protein_id'], s['protein_description'], s['protein_id_fasta'], s['AA_seq']
        fout_proteins_all.write('>{}\tunchanged\n{}\n'.format(protein_description, str(AA_seq)))
    
    conn.close()
    fout_proteins_changed.close()
    fout_proteins_all.close()
    
    # clear temp folder
    if keep_all:
        print('keep all intermediate files')
    else:
        print('remove temp folder')
        shutil.rmtree(tempfolder)

    print('finished!')

def runPerChomSqlite_vcf(file_mutations, file_sqlite, threads, outprefix, datatype, protein_keyword, keep_all, individual, chromosome_only, filter_PASS, chromosomes_genome, chromosomes_genome_description, file_gtf, info_field = None, info_field_thres=None):
    '''
    '''
    from vcf2mutation import convertVCF2MutationComplex
    from PrecisionProDB_vcf import readProtein2DF, openFile
    # get two mutation files from vcf file
    print('start extracting mutation file from the vcf input')
    outprefix_vcf = outprefix + '.vcf2mutation'
    individual = convertVCF2MutationComplex(file_vcf = file_mutations, outprefix = outprefix_vcf, individual_input=individual, filter_PASS = filter_PASS, chromosome_only = chromosome_only, info_field = info_field, info_field_thres=info_field_thres, threads = threads)
    individual = ','.join(individual)
    print('finished extracting mutations from the vcf file')
    file_mutations = outprefix + '.vcf2mutation.tsv'
    
    if TEST:
        print([file_sqlite, file_mutations, threads, outprefix, protein_keyword, datatype, keep_all, individual, chromosomes_genome, chromosomes_genome_description, file_gtf])
    
    runPerChomSqlite(file_sqlite, 
                     file_mutations, 
                     threads, 
                     outprefix, 
                     protein_keyword, 
                     datatype, 
                     keep_all, 
                     individual, 
                     chromosomes_genome, 
                     chromosomes_genome_description, 
                     file_gtf
                     )
    

    print('perGeno_vcf finished!')



def check_sqlite_file(file_path):
    """
    Checks if an SQLite file meets certain conditions.

    This function verifies the following conditions for the provided SQLite file:
    1. The file must exist.
    2. The file size must be greater than or equal to 100 bytes.
    3. The database must contain the tables 'CDSloc', 'protein_description', and 'chromosomes_using'.
    4. The database must contain at least one table name that starts with 'genomicLocs'.

    Args:
        file_path (str): The path to the SQLite file.

    Returns:
        bool: True if all conditions are met, otherwise False.
    """
    # Check if the file exists
    if not os.path.exists(file_path):
        return False
    
    # Check if the file size is greater than or equal to 100 bytes
    if os.path.getsize(file_path) < 100:
        return False
    
    try:
        # Connect to the SQLite database
        conn = sqlite3.connect(file_path)
        cursor = conn.cursor()

        # Get all table names from the database
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        table_names = [row[0] for row in cursor.fetchall()]
        
        # Check if required tables exist in the database
        required_tables = {"CDSloc", "protein_description", "chromosomes_using"}
        if not required_tables.issubset(table_names):
            return False

        # Check if at least one table name starts with "genomicLocs"
        if not any(table_name.startswith("genomicLocs") for table_name in table_names):
            return False

        return True
    
    except sqlite3.Error:
        # If there's an error with connecting to or querying the SQLite database, return False
        return False
    
    finally:
        # Close the database connection if it was opened
        if 'conn' in locals():
            conn.close()
            
def main_PrecsionProDB_Sqlite(file_genome, file_gtf, file_mutations, file_protein, threads, outprefix, datatype, protein_keyword, filter_PASS, individual, chromosome_only, keep_all, file_sqlite, info_field=None, info_field_thres=None):
    '''
    the main function of PrecisionProDB_Sqlite
    '''

    if os.path.exists(file_sqlite):
        if not check_sqlite_file(file_sqlite):
            print(f'sqlite file "{file_sqlite}" is not good? delete the file before running the program!!!')
            sys.exit()
        else:
            print('running in sqlite mode')
            print(f'use existing sqlite file "{file_sqlite}" for gene annotation')
    else:
        print('running in sqlite mode')
        print(f'sqlite file "{file_sqlite}" does not exist, create one first')
        buildSqlite.create_sqlite(file_sqlite, file_genome, file_gtf, file_protein, outprefix, datatype, protein_keyword, threads, keep_all)

    tempfolder = outprefix +'_temp'
    if not os.path.exists(tempfolder):
        os.mkdir(tempfolder)
    # get table chromosomes_using from file_sqlite as a dataframe table
    conn = buildSqlite.get_connection(file_sqlite)
    query = 'SELECT * FROM chromosomes_using'
    df_chromosomes = pd.read_sql_query(query, conn)
    conn.close()
    chromosomes_genome_description = list(df_chromosomes['chromosomes_genome_description'])
    chromosomes_genome = list(df_chromosomes['chromosome'])

    # run if file_mutations is a single variant
    pattern = re.compile(r'(chr)?(\d+)-(\d+)-([A-Za-z]+)-([A-Za-z]+)')
    match = pattern.match(file_mutations)
    if match and (not os.path.exists(file_mutations)):
        df_mutations = perChrom.parse_mutation(file_mutations=file_mutations)
        ls_results = []
        for chr_temp, file_mutations_single_chromosome in df_mutations.groupby('chr'):
            # chr_temp = df_mutations.iloc[0]['chr']
            chr_temp = get_k_new(chr_temp, chromosomes_genome, chromosomes_genome_description)
            if chr_temp in chromosomes_genome:
                chromosome = chr_temp
                file_mutations_single_chromosome['chr'] = chromosome
            else:
                print('chromosome not found in the genome')
                continue
            perchrom_sqlite = perChromSqlite.PerChrom_sqlite(file_sqlite = file_sqlite,
                            file_mutations = file_mutations_single_chromosome,
                            threads = threads,
                            outprefix = outprefix,
                            datatype = datatype,
                            chromosome = chromosome)
            print('run perChrom for chromosome', chromosome)
            ls_results.append(perchrom_sqlite.run_perChrom(save_results=False))
        df_transcript3 = pd.concat(ls_results, ignore_index=True)
        perChrom.save_mutation_and_proteins(df_transcript3, outprefix)
        # clear temp folder
        if keep_all:
            print('keep all intermediate files')
        else:
            print('remove temp folder')
            shutil.rmtree(tempfolder)
        
        # rename files
        if os.path.exists(outprefix + '.aa_mutations.csv'):
            os.rename(outprefix + '.aa_mutations.csv',outprefix + '.pergeno.aa_mutations.csv')
        if os.path.exists(outprefix + '.mutated_protein.fa'):
            os.rename(outprefix + '.mutated_protein.fa',outprefix + '.pergeno.mutated_protein.fa')
        print('finished!')
        
    elif file_mutations.lower().endswith('.vcf') or file_mutations.lower().endswith('.vcf.gz'):
        print('variant file is a vcf file')
        runPerChomSqlite_vcf(
            file_mutations = file_mutations, 
            file_sqlite = file_sqlite,
            threads=threads, 
            outprefix=outprefix, 
            datatype=datatype, 
            protein_keyword=protein_keyword, 
            keep_all=keep_all,
            individual=individual,
            filter_PASS = filter_PASS, 
            chromosome_only = chromosome_only,
            chromosomes_genome=chromosomes_genome,
            chromosomes_genome_description=chromosomes_genome_description,
            file_gtf=file_gtf,
            info_field=info_field,
            info_field_thres=info_field_thres
            )
    else:
        print('variant file is a tsv file')
        runPerChomSqlite(
            file_mutations = file_mutations, 
            file_sqlite = file_sqlite,
            threads=threads, 
            outprefix=outprefix, 
            datatype=datatype, 
            protein_keyword=protein_keyword, 
            keep_all=keep_all,
            individual=individual,
            chromosomes_genome=chromosomes_genome,
            chromosomes_genome_description=chromosomes_genome_description,
            file_gtf=file_gtf
            )



description = '''PrecisionProDB_Sqlite, personal proteogenomic tools which outputs a new reference protein based on the variants data. Use sqlite3 to store the gene annotations.
if datatype is gtf or not set, the gtf input is required
'''

def main():
    import argparse
    parser = argparse.ArgumentParser(description=description)
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-g','--genome', help = 'the reference genome sequence in fasta format. It can be a gzip file', default='')
    parser.add_argument('-f', '--gtf', help='gtf file with CDS and exon annotations. It can be a gzip file', default='')
    parser.add_argument('-m', '--mutations', help='''
                        a file stores the variants. 
                        If the file ends with ".vcf" or ".vcf.gz", treat as vcf input. Otherwise, treat as TSV input. 
                        a string like "chr1-788418-CAG-C" or "chr1-942451-T-C,1-6253878-C-T,1-2194700-C-G" can used as variant input, too. In this mode, --sample will not be used.
                        If multiple vcf files are provided, use "," to join the file names. For example, "--mutations file1.vcf,file2.vcf". A pattern match is also supported for input vcf, but quote is required to get it work. For example '--mutations "file*.vcf" ' 
                        
                        ''', default = '', required=False)
    parser.add_argument('-p','--protein', help = 'protein sequences in fasta format. It can be a gzip file. Only proteins in this file will be checked', default='')
    parser.add_argument('-t', '--threads', help='number of threads/CPUs to run the program. default, use all CPUs available', type=int, default=min(20, os.cpu_count()))
    parser.add_argument('-o', '--out', help='''output prefix, folder path could be included. Three or five files will be saved depending on the variant file format. Outputs include the annotation for mutated transcripts, the mutated or all protein sequences, two variant files from vcf. {out}.pergeno.aa_mutations.csv, {out}.pergeno.protein_all.fa, {out}.protein_changed.fa, {out}.vcf2mutation_1/2.tsv. default "perGeno" ''', default="perGeno")
    parser.add_argument('-a', '--datatype', help='''input datatype, could be GENCODE_GTF, GENCODE_GFF3, RefSeq, Ensembl_GTF or gtf. default "gtf". Ensembl_GFF3 is not supported. ''', default='gtf', type=str, choices=['GENCODE_GTF', 'GENCODE_GFF3','RefSeq','Ensembl_GTF','gtf'])
    parser.add_argument('-k','--protein_keyword', help='''field name in attribute column of gtf file to determine ids for proteins. default "auto", determine the protein_keyword based on datatype. "transcript_id" for GENCODE_GTF, "protein_id" for "RefSeq" and "Parent" for gtf and GENCODE_GFF3 ''', default='auto')
    parser.add_argument('-F', '--no_filter', help='default only keep variant with value "PASS" FILTER column of vcf file. if set, do not filter', action='store_true')
    parser.add_argument('-s', '--sample', help='''
                        sample name in the vcf to extract the variant information. default: None, extract the first sample. 
                        For multiple samples, use "," to join the sample names. For example, "--sample sample1,sample2,sample3". 
                        To use all samples, use "--sample ALL_SAMPLES". 
                        To use all variants regardless where the variants from, use "--sample ALL_VARIANTS".
                        ''', default=None)
    parser.add_argument('-A','--all_chromosomes', help='default keep variant in chromosomes and ignore those in short fragments of the genome. if set, use all chromosomes including fragments when parsing the vcf file', action='store_true')
    parser.add_argument('--keep_all', help='If set, do not delete files generated during the run', action='store_true')
    parser.add_argument('-S','--sqlite', help='''A path of sqlite file for re-use of annotation info. default '', do not use sqlite. The program will create a sqlite file if the file does not exist. If the file already exists, the program will use data stored in the file. It will cause error if the content in the sqlite file is not as expected. ''', default='', type=str)
    parser.add_argument('--info_field', help='fields to use in the INFO column of the vcf file to filter variants. Default None', default = None)
    parser.add_argument('--info_field_thres', help='threhold for the info field. Default None, do not filter any variants. If set "--info_filed AF --info_field_thres 0.01", only keep variants with AF >= 0.01', default = None)

    
    if TEST:
        f = parser.parse_args(f"-g {file_genome} -f {file_gtf} -m {file_mutations} -p {file_protein} -t {threads} -o {outprefix} -a {datatype} -k {protein_keyword} -F --keep_all -S {file_sqlite}".split())
    else:
        f = parser.parse_args()
    
    file_genome = f.genome
    file_gtf = f.gtf
    file_mutations = f.mutations
    file_protein = f.protein
    threads = f.threads
    outprefix = f.out
    datatype = f.datatype
    protein_keyword = f.protein_keyword
    filter_PASS = not f.no_filter
    individual = f.sample
    chromosome_only = not f.all_chromosomes
    keep_all = f.keep_all
    file_sqlite = f.sqlite
    print(description)
    print(f)

    main_PrecsionProDB_Sqlite(file_genome, file_gtf, file_mutations, file_protein, threads, outprefix, datatype, protein_keyword, filter_PASS, individual, chromosome_only, keep_all, file_sqlite, info_field=f.info_field, info_field_thres=f.info_field_thres)

if __name__ == '__main__':
    main()