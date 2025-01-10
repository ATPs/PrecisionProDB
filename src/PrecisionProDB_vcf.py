import pandas as pd
from Bio import SeqIO
import gzip
import os
from PrecisionProDB_core import PerGeno
from vcf2mutation import getMutationsFromVCF

def openFile(filename):
    '''open text or gzipped text file
    '''
    if filename.endswith('.gz'):
        return gzip.open(filename,'rt')
    return open(filename)


def readProtein2DF(filename):
    '''given a fasta file, return a dataframe with seq_id, seq_anno, and seq
    '''
    ls = list(SeqIO.parse(openFile(filename),'fasta'))
    if len(ls) == 0:
        return pd.DataFrame()
    tdf = pd.DataFrame()
    tdf['seq_id'] = [e.id for e in ls]
    tdf['seq_anno'] = [e.description for e in ls]
    tdf['seq_anno'] = tdf.apply(lambda x:x['seq_anno'].split(maxsplit=1)[-1] if x['seq_anno'].split(maxsplit=1)[-1] != x['seq_id'] else '', axis=1)
    tdf['seq'] = [str(e.seq) for e in ls]
    return tdf

def runPerGenoVCF(
                    file_genome,
                    file_gtf,
                    file_mutations,
                    file_protein,
                    threads = os.cpu_count(),
                    outprefix = 'perGeno',
                    datatype = 'gtf',
                    protein_keyword = 'auto',
                    filter_PASS = True,
                    individual = None,
                    chromosome_only = True,
                    keep_all = False
                ):
    '''
    run perGeno with a vcf file as variant input
    '''
    # get two mutation files from vcf file
    print('start extracting mutation file from the vcf input')
    outprefix_vcf = outprefix + '.vcf2mutation'
    getMutationsFromVCF(file_vcf = file_mutations, outprefix = outprefix_vcf, individual=individual, filter_PASS = filter_PASS, chromosome_only = chromosome_only)
    print('finished extracting mutations from the vcf file')
    file_mutations_1 = outprefix_vcf + '_1.tsv'
    file_mutations_2 = outprefix_vcf + '_2.tsv'

    # run perGeno for mutations_1
    print('start running PrecisionProDB for first strand of the genome mutation file')
    outprefix_1 = outprefix + '_1'
    pergeno_1 = PerGeno(file_genome = file_genome, file_gtf=file_gtf, file_mutations = file_mutations_1, file_protein=file_protein, threads=threads, outprefix=outprefix_1, datatype=datatype, protein_keyword=protein_keyword, keep_all=keep_all)
    print(pergeno_1.__dict__)
    pergeno_1.splitInputByChromosomes()
    pergeno_1.runPerChom()

    # run perGeno for mutations_2
    print('start running PrecisionProDB for second strand of the genome mutation file')
    outprefix_2 = outprefix + '_2'
    pergeno_2 = PerGeno(file_genome = file_genome, file_gtf=file_gtf, file_mutations = file_mutations_2, file_protein=file_protein, threads=threads, outprefix=outprefix_2, datatype=datatype, protein_keyword=protein_keyword, keep_all=keep_all)
    print(pergeno_1.__dict__)
    pergeno_2.splitInputByChromosomes()
    pergeno_2.runPerChom()

    # read in the changed protein sequences
    fout_protein_changed_1 = outprefix_1 + '.pergeno.protein_changed.fa'
    fout_protein_changed_2 = outprefix_2 + '.pergeno.protein_changed.fa'
    fout_protein_changed = outprefix +  '.pergeno.protein_changed.fa'
    df_seqs_changed_1 = readProtein2DF(fout_protein_changed_1)
    df_seqs_changed_2 = readProtein2DF(fout_protein_changed_2)
    df_seqs_changed_1['batch'] = '1'
    df_seqs_changed_2['batch'] = '2'
    df_seqs_changed = pd.concat([df_seqs_changed_1,df_seqs_changed_2], ignore_index=True)

    # read in the mutation files
    fout_mutations_1 = outprefix_1 + '.pergeno.aa_mutations.csv'
    fout_mutations_2 = outprefix_2 + '.pergeno.aa_mutations.csv'
    fout_mutations = outprefix + '.pergeno.aa_mutations.csv'
    try:
        df_mutations_1 = pd.read_csv(fout_mutations_1, sep='\t')
    except:
        df_mutations_1 = pd.DataFrame()
    try:
        df_mutations_2 = pd.read_csv(fout_mutations_2, sep='\t')
    except:
        df_mutations_2 = pd.DataFrame()
    df_mutations_1['batch'] = '1'
    df_mutations_2['batch'] = '2'
    df_mutations = pd.concat([df_mutations_1, df_mutations_2], ignore_index=True)

    # join the mutation and protein file
    df_joined = df_seqs_changed.merge(df_mutations, left_on = ['seq_id','batch'], right_on = ['protein_id_fasta','batch'])

    df_joined_dedup = df_joined.fillna('')
    df_joined_dedup = df_joined_dedup.groupby(by = [e for e in df_joined_dedup.columns if e != 'batch'])['batch'].apply(lambda x:';'.join(str(e) for e in x)).reset_index()

    print('total number of changed proteins from two strands of chromosome: {}. After removing protein seq_id with the same mutations, {} sequences left. After remove protein seq_id with same final sequence, {} sequences left.'.format(df_seqs_changed_1.shape[0] + df_seqs_changed_2.shape[0], df_joined_dedup.shape[0], df_joined.groupby(['seq_id','seq']).ngroups))

    protein_ids_changed = set(df_joined_dedup['protein_id_fasta'])
    protein_ids_changedBothStrands = set(df_seqs_changed_1['seq_id']) & set(df_seqs_changed_2['seq_id'])
    n_changed_both = len(protein_ids_changedBothStrands)
    n_changed_bothDiffMut = df_joined_dedup.shape[0] - len(protein_ids_changed)
    print('proteins changed in at least one strand: {}. Proteins changed in both strands: {}. Proteins changed in both strands and the mutation in both strands are identical: {}. Proteins changed in both strands, the final protein sequences are identical but the mutations are different: {}'.format(len(protein_ids_changed), n_changed_both, n_changed_both - n_changed_bothDiffMut, df_joined_dedup.shape[0] - df_joined.groupby(['seq_id','seq']).ngroups))
    
    ## save the result
    df_joined_dedup['seq_id_final'] = df_joined_dedup.apply(lambda r:r['seq_id'] + '__'+ r['batch'] if r['batch'] in ['1','2'] else r['seq_id'] +'__12',axis=1)
    f = open(fout_protein_changed,'w')
    for _, r in df_joined_dedup.iterrows():
        f.write('>{}\t{}\n{}\n'.format(r['seq_id_final'], r['seq_anno'], r['seq']))
    f.close()

    # combine all protein_sequences
    fout_protein_all_1 = outprefix_1 + '.pergeno.protein_all.fa'
    fout_protein_all_2 = outprefix_2 + '.pergeno.protein_all.fa'
    fout_protein_all = outprefix + '.pergeno.protein_all.fa'
    f = open(fout_protein_all,'w')
    for s in SeqIO.parse(openFile(file_protein),'fasta'):
        if s.id not in protein_ids_changedBothStrands:
            f.write('>{}\n{}\n'.format(s.description + '\tunchanged', str(s.seq)))
    f.write(open(fout_protein_changed).read())
    f.close()


    df_mutations = df_joined_dedup[df_mutations_1.columns].copy()
    df_mutations['protein_id_fasta'] = df_joined_dedup['seq_id_final']
    df_mutations.to_csv(fout_mutations, sep='\t',index=None)

    # clean up files
    if not keep_all:
        os.remove(fout_mutations_1)
        os.remove(fout_mutations_2)
        os.remove(fout_protein_changed_1)
        os.remove(fout_protein_changed_2)
        os.remove(fout_protein_all_1)
        os.remove(fout_protein_all_2)
    print('perGeno_vcf finished!')


description = '''PrecisionProDB_vcf, personal proteogenomic tools which outputs a new reference protein based on the variants data. VCF file as the variant input
'''

def main():
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-g','--genome', help = 'the reference genome sequence in fasta format. It can be a gzip file', required=True)
    parser.add_argument('-f', '--gtf', help='gtf file with CDS and exon annotations. It can be a gzip file', required=True)
    parser.add_argument('-m', '--mutations', help='a file stores the variants', required=True)
    parser.add_argument('-p','--protein', help = 'protein sequences in fasta format. It can be a gzip file. Only proteins in this file will be checked', required=True)
    parser.add_argument('-t', '--threads', help='number of threads/CPUs to run the program. default, use all CPUs available', type=int, default=min(20, os.cpu_count()))
    parser.add_argument('-o', '--out', help='''output prefix. Five files will be saved, including the annotation for mutated transcripts, the mutated or all protein sequences, two variant files from vcf. {out}.pergeno.aa_mutations.csv, {out}.pergeno.protein_all.fa, {out}.protein_changed.fa, {out}.vcf2mutation_1/2.tsv. default "perGeno" ''', default="perGeno")
    parser.add_argument('-a', '--datatype', help='''input datatype, could be GENCODE_GTF, GENCODE_GFF3, RefSeq, Ensembl_GTF or gtf. default "gtf". Ensembl_GFF3 is not supported. ''', default='gtf', type=str, choices=['GENCODE_GTF', 'GENCODE_GFF3','RefSeq','Ensembl_GTF','gtf'])
    parser.add_argument('-k','--protein_keyword', help='''field name in attribute column of gtf file to determine ids for proteins. default "auto", determine the protein_keyword based on datatype. "transcript_id" for GENCODE_GTF, "protein_id" for "RefSeq" and "Parent" for gtf and GENCODE_GFF3 ''', default='auto')
    parser.add_argument('-F', '--no_filter', help='default only keep variant with value "PASS" FILTER column of vcf file. if set, do not filter', action='store_true')
    parser.add_argument('-s', '--sample', help='sample name in the vcf to extract the variant information. default: None, extract the first sample', default=None)
    parser.add_argument('-A','--all_chromosomes', help='default keep variant in chromosomes and ignore those in short fragments of the genome. if set, use all chromosomes including fragments when parsing the vcf file', action='store_true')
    parser.add_argument('--keep_all', help='If set, do not delete files generated during the run', action='store_true')

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

    print(f)

    runPerGenoVCF(
        file_genome = file_genome, 
        file_gtf=file_gtf, 
        file_mutations = file_mutations, 
        file_protein=file_protein, 
        threads=threads, 
        outprefix=outprefix, 
        datatype=datatype, 
        protein_keyword=protein_keyword, 
        filter_PASS=filter_PASS, 
        individual=individual, 
        chromosome_only=chromosome_only, 
        keep_all=keep_all
        )

if __name__ == '__main__':
    main()