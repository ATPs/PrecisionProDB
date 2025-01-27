import gzip
import re
import io
import pandas as pd
import glob
import os
import itertools
from multiprocessing import Pool
import sqlite3
import numpy as np

try:
    import tqdm
except:
    print('Cannot import tqdm. Will not use tqdm.')
    tqdm = False

import itertools

def grouper_it(n, iterable, M=1):
    """
    Groups an iterable into chunks of size `n` and further groups these chunks into outer lists of size `M`.

    This function is memory-efficient and works lazily, making it suitable for large or infinite iterables.

    Parameters:
        n (int): The size of each inner chunk. Each chunk will contain up to `n` items.
        iterable (iterable): The input iterable (e.g., list, tuple, generator) to be grouped.
        M (int, optional): The number of chunks to group into each outer list. Defaults to 1.

    Yields:
        list of lists: Each yielded value is a list containing `M` chunks, where each chunk is a list of up to `n` items.
                       If the iterable is exhausted, the last yielded value may contain fewer than `M` chunks.

    Examples:
        >>> my_list = list(range(1, 21))  # [1, 2, 3, ..., 20]
        >>> for outer_list in grouper_it(3, my_list, M=2):
        ...     print(outer_list)
        [[1, 2, 3], [4, 5, 6]]
        [[7, 8, 9], [10, 11, 12]]
        [[13, 14, 15], [16, 17, 18]]
        [[19, 20]]

        >>> for outer_list in grouper_it(2, my_list, M=3):
        ...     print(outer_list)
        [[1, 2], [3, 4], [5, 6]]
        [[7, 8], [9, 10], [11, 12]]
        [[13, 14], [15, 16], [17, 18]]
        [[19, 20]]

        >>> for outer_list in grouper_it(4, my_list, M=1):
        ...     print(outer_list)
        [[1, 2, 3, 4]]
        [[5, 6, 7, 8]]
        [[9, 10, 11, 12]]
        [[13, 14, 15, 16]]
        [[17, 18, 19, 20]]

    Notes:
        - If `n` or `M` is zero, the function will yield empty lists or no lists, depending on the input.
        - If the iterable is exhausted before filling `M` chunks, the last yielded value will contain the remaining chunks.
        - This function is memory-efficient because it uses `itertools.islice` to process the iterable lazily.
    """
    it = iter(iterable)
    while True:
        # Collect M chunks, each of size n
        outer_chunk = []
        for _ in range(M):
            chunk = list(itertools.islice(it, n))
            if not chunk:  # Stop if no more items are left
                break
            outer_chunk.append(chunk)
        if not outer_chunk:  # Stop if no chunks were collected
            break
        yield outer_chunk


CHROMOSOMES = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'}

def openFile(filename):
    '''open text or gzipped text file
    '''
    if filename.endswith('.gz'):
        return gzip.open(filename,'rt')
    return open(filename)

def getMutationsFromVCF(file_vcf, outprefix = None, individual=None, filter_PASS = True, chromosome_only = True):
    '''return two dataframe of mutations.
    if outprefix is not None, write the two dataframe to outprfex +'_1/2.tsv' two files as mutation file
    file_vcf is the vcf file, it could be a gzip file
    if individual is None, work on the first individual in the vcf file
    default, only work on variants in chromosome
    '''
    if ',' in file_vcf:
        files_vcf = file_vcf.split(',')
    elif os.path.exists(file_vcf):
        files_vcf = [file_vcf]
    else:
        files_vcf = glob.glob(file_vcf)
        if len(files_vcf) == 0:
            print(file_vcf, 'not found.')
            return None
    
    ls_keep = []
    write_header = True
    
    for file_vcf in files_vcf:
        fo = openFile(file_vcf)
        for line in fo:
            if not line.startswith('##'):
                break
        
        columns = line.strip().split('\t')
        # get the column to keep in the 
        if individual is None:
            indi_col = 9
            print('individual is None, select the first sample in vcf file, which is', columns[9])
        else:
            if individual in columns:
                indi_col = columns.index(individual)
                print('individual is provided as', individual, 'which exists in the vcf file')
            else:
                print(individual, 'does not exist in the vcf file')
                exit(-1)
        
        # chromosomes
        chromosomes = [str(i) for i in range(1,23)] + list('XY')
        chromosomes = ['chr' + i for i in chromosomes] + chromosomes
        chromosomes = set(chromosomes)
        
        if write_header:
            ls_keep.append('chr\tpos\tref\talt1\talt2\n')
            write_header = False
        
        for line in fo:
            es = line.strip().split('\t')
            chromosome, position, reference, alternatives, filter, genotype = es[0], es[1], es[3], es[4], es[6], es[indi_col]

            if chromosome_only:
                # skip if not chromosome but some scaffolds
                if chromosome not in chromosomes:
                    continue
            
            if filter_PASS:
                # skip if not "PASS"
                if filter != "PASS":
                    continue

            # skip if no mutation
            if genotype.startswith('0|0') or genotype.startswith('0/0'):
                continue

            alternatives = alternatives.split(',')
            genotype = genotype.split(':')[0]
            if genotype == '.':
                continue
            
            if '|' in genotype:
                GTs = genotype.split('|')
            else:
                GTs = genotype.split('/')
            GTs = [int(e) if e != '.' else 0 for e in GTs]
            alternatives = [reference] + alternatives
            alternative1, alternative2 = alternatives[GTs[0]], alternatives[GTs[1]]
            ls_keep.append('{}\t{}\t{}\t{}\t{}\n'.format(chromosome, position, reference, alternative1, alternative2))
    
    df = pd.read_csv(io.StringIO(''.join(ls_keep)), sep='\t',low_memory=False)
    print('before QC, sites with mutations:', df.shape[0],'of which, homozyous sites:', df[df['alt1'] == df['alt2']].shape[0])

    df1 = df[['chr', 'pos', 'ref', 'alt1']].copy()
    df2 = df[['chr', 'pos', 'ref', 'alt2']].copy()
    column_names = ['chr', 'pos', 'ref', 'alt']
    df1.columns = column_names
    df2.columns = column_names

    # filter dataframe
    df1 = df1[(df1['ref'] != df1['alt']) & (df1['alt'] != '*')]
    df2 = df2[(df2['ref'] != df2['alt']) & (df2['alt'] != '*')]
    df1 = df1.drop_duplicates(subset=['chr','pos','ref'], keep='first')
    df2 = df2.drop_duplicates(subset=['chr','pos','ref'], keep='first')
    print('after QC, number of variants in two dataframes are:', df1.shape[0], df2.shape[0])

    # save results
    if outprefix is not None:
        fout1 = outprefix + '_1.tsv'
        fout2 = outprefix + '_2.tsv'
        df1.to_csv(fout1, sep='\t', index=None)
        df2.to_csv(fout2, sep='\t', index=None)
        return [fout1, fout2]

    return df1, df2


def split_complex_variant(ls_to_write):
    '''
    ls_to_write = [chromosome, position, reference, alternative, '\t'.join(str(e) for e in ls_alternatives)]
    or ls_to_write = [chromosome, position, reference, alternative]
    deal with complex variants. usually it should be like chr1-111-A-T, chr1-111-AT-A, chr1-111-A-AT. but some times, it could be like chr1-111-GG-TT
    '''
    chromosome, position, reference, alternative = ls_to_write[:4]
    position = int(position)
    ls = []
    n = min(len(reference), len(alternative))
    for i in range(n-1):
        ref = reference[i]
        alt = alternative[i]
        if ref != alt:
            ls.append([chromosome, position + i, ref, alt] + ls_to_write[4:])
    ref = reference[n-1:]
    alt = alternative[n-1:]
    ls.append([chromosome, position + n-1, ref, alt] + ls_to_write[4:])
    return ls

def get_writing_string_for_complex_variant(ls_to_write):
    '''
    ls_to_write = [chromosome, position, reference, alternative, '\t'.join(str(e) for e in ls_alternatives)]
    or ls_to_write = [chromosome, position, reference, alternative]
    deal with complex variants. usually it should be like chr1-111-A-T, chr1-111-AT-A, chr1-111-A-AT. but some times, it could be like chr1-111-GG-TT
    
    '''
    ls = split_complex_variant(ls_to_write)
    ls_txt = ['\t'.join([str(e) for e in i]) + '\n' for i in ls]
    return ''.join(ls_txt)


def vcfINFO2dict(info_str):
    """
    Parse a VCF INFO string into a dictionary, preserving all values as strings.
    
    Args:
        info_str (str): The INFO field from a VCF file (e.g., "DP=100;DB;AA=T")
        
    Returns:
        dict: Dictionary where keys are INFO tags and values are strings
              (empty string for flags without values)
    """
    info_dict = {}
    
    if not info_str or info_str == ".":
        return info_dict
        
    for entry in info_str.split(";"):
        if "=" in entry:
            key, value = entry.split("=", 1)
            try:
                value = float(value)
            except:
                pass
            info_dict[key] = value
        else:
            # Handle flag-type entries without values
            info_dict[entry] = 1
            
    return info_dict


def processOneLineOfVCF(line, 
                        individual_col, 
                        chromosome_only=True, 
                        filter_PASS=True, 
                        info_field=None, 
                        info_field_thres=None, 
                        chromosomes=CHROMOSOMES,
                        file_vcf=None
                        ):
    '''Process one line of VCF file and return formatted mutation information
    
    Args:
        line (str): One line from VCF file
        individual_col (list): List of column indices for individuals/samples
        chromosome_only (bool): If True, only process chromosomes (skip scaffolds)
        filter_PASS (bool): If True, only process variants with FILTER=PASS
        info_field (str): INFO field key to filter by (e.g. 'AF')
        info_field_thres (float): Threshold for info_field filtering
        chromosomes (set): Set of valid chromosome names
        file_vcf (str): VCF file path for error reporting
        
    Returns:
        str: Formatted mutation line(s) or empty string if variant should be skipped
    '''
    # print(line)
    es = line.strip().split('\t')
    chromosome, position, reference, alternatives, FILTER = es[0], es[1], es[3], es[4], es[6]
    genotypes = [es[i] for i in individual_col]

    if len(es) < 8:
        print(f'Warning: Malformed VCF line (too few columns) in {file_vcf}: {line}')

    if chromosome_only:
        # skip if not chromosome but some scaffolds
        if chromosome not in chromosomes:
            return ''
    
    if filter_PASS:
        # skip if not "PASS"
        if FILTER != "PASS":
            return ''

    # filter by info field
    if info_field:
        if info_field_thres:
            info_str = es[7]
            info_dict = vcfINFO2dict(info_str)
            if info_field in info_dict:
                if info_dict[info_field] < info_field_thres:
                    return ''
            else:
                print('waring! INFO field', info_field, 'not found in', file_vcf, line)
    # skip if no mutation
    if all([genotype.startswith('0|0') or genotype.startswith('0/0') or genotype.startswith('.|.') or genotype.startswith('./.') for genotype in genotypes]) and len(genotypes) >= 1:
        return ''

    alternatives = alternatives.split(',')
    genotypes = [genotype.split(':')[0] for genotype in genotypes]
    # if '.' in genotypes, change to './.'
    genotypes = ['./.' if genotype == '.' else genotype for genotype in genotypes]
    
    if all(['.' in genotype for genotype in genotypes]) and len(genotypes) >= 1:
        print('line with too many missing genotypes skipped')
        return ''
    
    GTs = [genotype.split('|') if '|' in genotype else genotype.split('/') for genotype in genotypes]
    GTs = [[int(e) if e != '.' else 0 for e in i] for i in GTs]
    
    alleles = [reference] + alternatives
    ls_new_line = []
    for alternative_index in range(len(alternatives)):
        alternative = alternatives[alternative_index]
        ls_alternatives = []
        for GT in GTs:
            for i in GT:
                if i == alternative_index + 1:
                    ls_alternatives.append(1)
                else:
                    ls_alternatives.append(0)
        if ls_alternatives:
            ls_to_write = [chromosome, position, reference, alternative, '\t'.join(str(e) for e in ls_alternatives)]
        else:
            ls_to_write = [chromosome, position, reference, alternative]
        new_line = get_writing_string_for_complex_variant(ls_to_write)
        ls_new_line.append(new_line)
    
    return ''.join(ls_new_line)

def processManyLineOfVCF(lines, 
                        individual_col, 
                        chromosome_only=True, 
                        filter_PASS=True, 
                        info_field=None, 
                        info_field_thres=None, 
                        chromosomes=CHROMOSOMES,
                        file_vcf=None
                        ):
    
    return ''.join(processOneLineOfVCF(line, individual_col, chromosome_only, filter_PASS, info_field, info_field_thres, chromosomes, file_vcf) for line in lines)



def sanitize_identifier(identifier: str) -> str:
    """
    Sanitize a string to make it a valid identifier (e.g., for SQLite table names, column names, etc.).

    Args:
        identifier (str): The original string to sanitize.

    Returns:
        str: A sanitized string suitable for use as an identifier.
    """
    # Replace invalid characters with underscores
    sanitized_name = re.sub(r'[^a-zA-Z0-9_]', '_', identifier)
    
    # Ensure the name starts with a letter or underscore
    if not sanitized_name[0].isalpha() and sanitized_name[0] != '_':
        sanitized_name = '_' + sanitized_name
    
    return sanitized_name


def tsv2sqlite(tsv_file: str, sqlite_file: str = None, batch_size: int = 100):
    """
    Convert a TSV file to a SQLite database, using the row index as the primary key.

    Args:
        tsv_file (str): Path to the input TSV file.
        sqlite_file (str): Path to the output SQLite database file.
        batch_size (int): Number of rows to process in each batch. Default is 100.
    """
    if sqlite_file is None:
        sqlite_file = tsv_file + '.sqlite'
    
    # remove sqlite file if it already exists
    if os.path.exists(sqlite_file):
        print(f"removing existing sqlite file: {sqlite_file}")
        os.remove(sqlite_file)
    
    # Create a connection to the SQLite database
    conn = sqlite3.connect(sqlite_file)
    cursor = conn.cursor()

    # Read the TSV file in chunks
    reader = pd.read_csv(tsv_file, sep='\t', chunksize=batch_size)

    # Process each chunk
    print('storeing tsv file to sqlite database...')
    if tqdm:
        to_iter_value = tqdm.tqdm(enumerate(reader))
    else:
        to_iter_value = enumerate(reader)
    for i, chunk in to_iter_value:
        columns = [i for i in chunk.columns]
        column_group_N = len(columns) // 1000 + 1
        # split columns to column_group_N parts
        ls_columns = [columns[i::column_group_N] for i in range(column_group_N)]
        for table_name, columns_use in enumerate(ls_columns):
            table_name = 'table_'+str(table_name)
            if i == 0:
                # Generate the CREATE TABLE query dynamically based on the columns
                create_table_query = f'''
                CREATE TABLE IF NOT EXISTS {table_name} (
                    {', '.join([f'{col} TEXT' for col in columns_use])}
                );
                '''
                cursor.execute(create_table_query)
                conn.commit()

            # Insert the chunk into the SQLite table
            chunk[columns_use].to_sql(table_name, conn, if_exists='append', index=False)

    # Close the database connection
    conn.close()

    print(f"TSV file '{tsv_file}' has been successfully converted to SQLite database '{sqlite_file}'.")

def tsv2memmap(tsv_file, individuals = None, memmap_file=None, batch_size=100):
    '''tsv file is a tsv file with columns: chr, pos, ref, sample1__1, sample1__2, sample2__1, sample2__2, ..., sampleN__1, sampleN__2 columns
    individuals is a list of individuals to be used. the column values must be 0 or 1
    if memmap_file is None, memmap_file = tsv_file + '.memmap'
    '''
    if memmap_file is None:
        memmap_file = tsv_file + '.memmap'
        
    if individuals is None:
        tdf = pd.read_csv(tsv_file, sep='\t', nrows=1)
        individuals = [i for i in tdf.columns if i not in ['chr', 'pos', 'ref', 'alt','pos_end']]
    elif ',' in individuals:
        individuals = individuals.split(',')
    elif isinstance(individuals, str):
        individuals = [individuals]
    else:
        individuals = individuals
    
    # Read the TSV file in chunks
    reader = pd.read_csv(tsv_file, sep='\t', chunksize=batch_size, usecols=individuals, dtype='int8')

    # Get total number of rows
    total_rows = sum(1 for _ in open(tsv_file)) - 1  # Subtract header row
    
    # Initialize memory-mapped file
    shape = (total_rows, len(individuals))
    mmap = np.memmap(memmap_file, dtype='int8', mode='w+', shape=shape)
    
    # Process each chunk
    print('storing tsv file to memmap...')
    if tqdm:
        to_iter_value = tqdm.tqdm(enumerate(reader), total = total_rows//batch_size + 1)
    else:
        to_iter_value = enumerate(reader)

    # Write data to memmap
    row_offset = 0
    for i, chunk in to_iter_value:
        chunk = chunk[individuals]
        chunk_size = chunk.shape[0]
        mmap[row_offset:row_offset + chunk_size] = chunk.values
        row_offset += chunk_size
        
    # Flush changes to disk
    mmap.flush()
    del mmap
    print(f"TSV file '{tsv_file}' has been successfully converted to memory-mapped file '{memmap_file}'")


def convertVCF2MutationComplex(file_vcf, outprefix = None, individual_input="ALL_SAMPLES", filter_PASS = True, chromosome_only = True, info_field = None, info_field_thres=None, threads = 1):
    '''convert vcf file to tsv file. with columns: chr, pos, ref, sample1__1, sample1__2, sample2__1, sample2__2, ..., sampleN__1, sampleN__2 columns
    If individual is None, use the first sample in the vcf file. 
    If individual == 'ALL_SAMPLES', use all samples in the vcf file.
    If individual == 'ALL_VARIANTS', ignore sample columns
    info_field is the INFO field key used to filter. In ProHap, the default setting is "AF"
    info_field_thres is the threshold of the INFO field key used to filter. In ProHap, the default setting is 0.01
    
    '''
    if ',' in file_vcf:
        files_vcf = file_vcf.split(',')
    elif os.path.exists(file_vcf):
        files_vcf = [file_vcf]
    else:
        files_vcf = glob.glob(file_vcf)
        if len(files_vcf) == 0:
            print(file_vcf, 'not found.')
            return []
    
    if outprefix is None:
        file_output = file_vcf + '.tsv'
    else:
        file_output = outprefix + '.tsv'
        
    if info_field:
        if info_field_thres:
            info_field_thres = float(info_field_thres)
        
    fout = open(file_output, 'w')

    write_header = True
    for file_vcf in files_vcf:
        print('start converting vcf file:', file_vcf)
        fo = openFile(file_vcf)
        for line in fo:
            if not line.startswith('##'):
                break
        
        columns = line.strip().split('\t')
        # get the column to keep in the 
        # chromosomes
        chromosomes = CHROMOSOMES
        
        if individual_input == 'ALL_SAMPLES':
            individual = columns[9:]
            individual_col = [columns.index(indi) for indi in individual]
        elif individual_input == 'ALL_VARIANTS':
            individual_col = []
            individual = []
        elif individual_input is None:
            individual_col = [9]
            individual = [columns[9]]
        else:
            individual = list(dict.fromkeys([i for i in individual_input.split(',') if i])) # remove duplicate and empty values
            if any([i not in columns for i in individual]):
                print('some samples were not found in vcf file:', [i for i in individual if i not in columns])
                return None
            individual_col = [columns.index(indi) for indi in individual]
        
        column_names = ['chr', 'pos', 'ref', 'alt']
        column_for_samples = []
        for indi in individual_col:
            column_for_samples.append(columns[indi] + '__1')
            column_for_samples.append(columns[indi] + '__2')
        
        column_names = column_names + column_for_samples
        
        if write_header:
            fout.write('\t'.join(column_names) + '\n')
            write_header = False
        
        # print(threads)
        # use single thread is thread set to 1 or file_vcf smaller than 10MB
        if threads == 1 or os.path.getsize(file_vcf) < 10000000:
            for line in fo:
                new_line = processOneLineOfVCF(line, individual_col, chromosome_only, filter_PASS, info_field, info_field_thres, chromosomes, file_vcf)
                fout.write(new_line)
        else:
            pool = Pool(threads)
            if tqdm:
                to_iter_value = tqdm.tqdm(grouper_it(1000, fo, threads))
            else:
                to_iter_value = grouper_it(1000, fo, threads)
            for ls_lines in to_iter_value:
                new_lines = pool.starmap(processManyLineOfVCF, [(lines, individual_col, chromosome_only, filter_PASS, info_field, info_field_thres, chromosomes, file_vcf) for lines in ls_lines])
                fout.write(''.join(new_lines))
            pool.close()
            pool.join()
        
        print('finished converting vcf file:', file_vcf)
    fout.close()
    return column_for_samples

description = '''convert extract mutation information from vcf file
If --sample is not set or only one sample were used, the same as PrecisionProDB version 1.0. two files will be generated.
Otherwise, PrecisionProDB version 2.0 mode. Output a single file with extra columns for each sample. If only one sample is used and version 2.0 mode is desired, use "--sample sample,"
'''

def main():
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--file_vcf', help='''input vcf file. It can be a gzip file. If multiple vcf files are provided, use "," to join the file names. For example, "--file_vcf file1.vcf,file2.vcf". A pattern match is also supported, but quote is required to get it work. For example '--file_vcf "file*.vcf"' ''', required=True)
    parser.add_argument('-o', '--outprefix', help='output prefix to store the two output dataframes, default: None, do not write the result to files. file will be outprefix_1/2.tsv', default=None)
    parser.add_argument('-s', '--sample', help='sample name in the vcf to extract the variant information. default: None, extract the first sample. For multiple samples, use "," to join the sample names. For example, "--sample sample1,sample2,sample3". To use all samples, use "--sample ALL_SAMPLES". To use all variants regardless where the variants from, use "--sample ALL_VARIANTS".', default=None)
    parser.add_argument('-F', '--no_filter', help='default only keep variant with value "PASS" FILTER column of vcf file. if set, do not filter', action='store_true')
    parser.add_argument('-A','--all_chromosomes', help='default keep variant in chromosomes and ignore those in short fragments of the genome. if set, use all chromosomes including fragments', action='store_true')
    parser.add_argument('--info_field', help='fields to use in the INFO column of the vcf file to filter variants. Default None', default = None)
    parser.add_argument('--info_field_thres', help='threhold for the info field. Default None, do not filter any variants. If set "--info_filed AF --info_field_thres 0.01", only keep variants with AF >= 0.01', default = None)
    parser.add_argument('-t', '--threads', help='number of threads/CPUs to run the program. default, 1', type=int, default=1)

    f = parser.parse_args()
    chromosome_only = not f.all_chromosomes
    filter_PASS = not f.no_filter
    sample = f.sample
    #print(f)
    
    if sample is None:
        print('convert vcf to mutation information in version 1.0 mode')
        getMutationsFromVCF(file_vcf = f.file_vcf, outprefix = f.outprefix, individual=f.sample, filter_PASS = filter_PASS, chromosome_only = chromosome_only)
    elif ',' in  sample or sample == 'ALL_SAMPLES' or sample == 'ALL_VARIANTS':
        print('convert vcf to mutation information in version 2.0 mode')
        convertVCF2MutationComplex(file_vcf = f.file_vcf, outprefix = f.outprefix, individual_input=f.sample, filter_PASS = filter_PASS, chromosome_only = chromosome_only, info_field = f.info_field, info_field_thres = f.info_field_thres, threads = f.threads)
    else:
        print('convert vcf to mutation information in version 1.0 mode')
        getMutationsFromVCF(file_vcf = f.file_vcf, outprefix = f.outprefix, individual=f.sample, filter_PASS = filter_PASS, chromosome_only = chromosome_only)


if __name__ == '__main__':
    main()