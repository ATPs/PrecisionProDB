import gzip
import re
import io
import pandas as pd


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
    
    ls_keep = []
    ls_keep.append('chr\tpos\tref\talt1\talt2\n')
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


description = '''convert extract mutation information from vcf file
'''

def main():
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--file_vcf', help='input vcf file. It can be a gzip file', required=True)
    parser.add_argument('-o', '--outprefix', help='output prefix to store the two output dataframes, default: None, do not write the result to files. file will be outprefix_1/2.tsv', default=None)
    parser.add_argument('-s', '--sample', help='sample name in the vcf to extract the variant information. default: None, extract the first sample', default=None)
    parser.add_argument('-F', '--no_filter', help='default only keep variant with value "PASS" FILTER column of vcf file. if set, do not filter', action='store_true')
    parser.add_argument('-A','--all_chromosomes', help='default keep variant in chromosomes and ignore those in short fragments of the genome. if set, use all chromosomes including fragments', action='store_true')
    
    f = parser.parse_args()
    chromosome_only = not f.all_chromosomes
    filter_PASS = not f.no_filter
    #print(f)

    getMutationsFromVCF(file_vcf = f.file_vcf, outprefix = f.outprefix, individual=f.sample, filter_PASS = filter_PASS, chromosome_only = chromosome_only)


if __name__ == '__main__':
    main()