import gzip
import re
import io
import pandas as pd
import glob
import os


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


def convertVCF2MutationComlex(file_vcf, outprefix = None, individual="ALL_SAMPLES", filter_PASS = True, chromosome_only = True):
    '''convert vcf file to tsv file. with columns: chr, pos, ref, sample1__1, sample1__2, sample2__1, sample2__2, ..., sampleN__1, sampleN__2 columns
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
    
    if outprefix is None:
        file_output = file_vcf + '.tsv'
    else:
        file_output = outprefix + '.tsv'
    
    individual_input = individual
    
    fout = open(file_output, 'w')

    write_header = True
    for file_vcf in files_vcf:
        fo = openFile(file_vcf)
        for line in fo:
            if not line.startswith('##'):
                break
        
        columns = line.strip().split('\t')
        # get the column to keep in the 
        # chromosomes
        chromosomes = [str(i) for i in range(1,23)] + list('XY')
        chromosomes = ['chr' + i for i in chromosomes] + chromosomes
        chromosomes = set(chromosomes)
        
        
        if individual_input == 'ALL_SAMPLES':
            individual = columns[9:]
            individual_col = [columns.index(indi) for indi in individual]
        elif individual_input == 'ALL_VARIANTS':
            individual_col = []
            individual = []
        else:
            individual = list(dict.fromkeys([i for i in individual_input.split(',') if i])) # remove duplicate and empty values
            if any([i not in columns for i in individual]):
                print('some samples were not found in vcf file:', [i for i in individual if i not in columns])
                return None
            individual_col = [columns.index(indi) for indi in individual]
        
        column_names = ['chr', 'pos', 'ref', 'alt']
        for indi in individual_col:
            column_names.append(columns[indi] + '__1')
            column_names.append(columns[indi] + '__2')
        
        if write_header:
            fout.write('\t'.join(column_names) + '\n')
            write_header = False
        
        for line in fo:
            # print(line)
            es = line.strip().split('\t')
            chromosome, position, reference, alternatives, FILTER = es[0], es[1], es[3], es[4], es[6]
            genotypes = [es[i] for i in individual_col]

            if chromosome_only:
                # skip if not chromosome but some scaffolds
                if chromosome not in chromosomes:
                    continue
            
            if filter_PASS:
                # skip if not "PASS"
                if FILTER != "PASS":
                    continue

            # skip if no mutation
            if all([genotype.startswith('0|0') or genotype.startswith('0/0') or genotype.startswith('.|.') or genotype.startswith('./.') for genotype in genotypes]) and len(genotypes) >= 1:
                continue

            alternatives = alternatives.split(',')
            genotypes = [genotype.split(':')[0] for genotype in genotypes]
            # if '.' in genotypes, change to './.'
            genotypes = ['./.' if genotype == '.' else genotype for genotype in genotypes]
            
            if all(['.' in genotype for genotype in genotypes]) and len(genotypes) >= 1:
                print('line with too many missing genotypes skipped')
                continue
            
            GTs = [genotype.split('|') if '|' in genotype else genotype.split('/') for genotype in genotypes]
            GTs = [[int(e) if e != '.' else 0 for e in i] for i in GTs]
            
            alleles = [reference] + alternatives
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
                new_line = '\t'.join(ls_to_write) + '\n'
                fout.write(new_line)
        
    fout.close()
    return None

description = '''convert extract mutation information from vcf file
If --sample is not set or only one sample were used, the same as PrecisionProDB version 1.0. two files will be generated.
Otherwise, PrecisionProDB version 2.0 mode. Output a single file with extra columns for each sample. If only one sample is used and version 2.0 mode is desired, use "--sample sample,"
'''

def main():
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--file_vcf', help='''input vcf file. It can be a gzip file. If multiple vcf files are provided, use "," to join the file names. For example, "--file_vcf file1.vcf,file2.vcf" A pattern match is also supported, but quote is required to get it work. For example '--file_vcf "file*.vcf"' ''', required=True)
    parser.add_argument('-o', '--outprefix', help='output prefix to store the two output dataframes, default: None, do not write the result to files. file will be outprefix_1/2.tsv', default=None)
    parser.add_argument('-s', '--sample', help='sample name in the vcf to extract the variant information. default: None, extract the first sample. For multiple samples, use "," to join the sample names. For example, "--sample sample1,sample2,sample3". To use all samples, use "--sample ALL_SAMPLES". To use all variants regardless where the variants from, use "--sample ALL_VARIANTS".', default=None)
    parser.add_argument('-F', '--no_filter', help='default only keep variant with value "PASS" FILTER column of vcf file. if set, do not filter', action='store_true')
    parser.add_argument('-A','--all_chromosomes', help='default keep variant in chromosomes and ignore those in short fragments of the genome. if set, use all chromosomes including fragments', action='store_true')

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
        convertVCF2MutationComlex(file_vcf = f.file_vcf, outprefix = f.outprefix, individual=f.sample, filter_PASS = filter_PASS, chromosome_only = chromosome_only)
    else:
        print('convert vcf to mutation information in version 1.0 mode')
        getMutationsFromVCF(file_vcf = f.file_vcf, outprefix = f.outprefix, individual=f.sample, filter_PASS = filter_PASS, chromosome_only = chromosome_only)


if __name__ == '__main__':
    main()