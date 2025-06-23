from PrecisionProDB_core import PerGeno
from PrecisionProDB_vcf import runPerGenoVCF
import os
import downloadHuman
import time
import sys
import re
import buildSqlite

def get_version():
    """Read version from version file"""
    try:
        version_file = os.path.join(os.path.dirname(__file__), 'version')
        with open(version_file, 'r') as f:
            return f.read().strip()
    except FileNotFoundError:
        return "unknown"

description = '''
PrecisionProDB, a personal proteogenomic tool which outputs a new reference protein based on the variants data. 
A VCF or /a tsv file can be used as the variant input. If the variant file is in tsv format, at least four columns are required in the header: chr, pos, ref, alt. Additional columns will be ignored. Try to Convert the file to proper format if you have a bed file or other types of variant file. The pos column is 1-based like in the vcf file.
Additionally, a string like "chr1-788418-CAG-C" can used as variant input. It has to be combined with the --sqlite for quick check of the mutation effects
'''
def main():
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-v', '--version', action='version', version=f'PrecisionProDB {get_version()}')
    parser.add_argument('-g','--genome', help = 'the reference genome sequence in fasta format. It can be a gzip file', default='')
    parser.add_argument('-f', '--gtf', help='gtf file with CDS and exon annotations. It can be a gzip file', default='')
    parser.add_argument('-m', '--mutations', help='''
                        a file stores the variants. 
                        If the file ends with ".vcf" or ".vcf.gz", treat as vcf input. Otherwise, treat as TSV input. 
                        A string like "chr1-788418-CAG-C" or "chr1-942451-T-C,1-6253878-C-T,1-2194700-C-G" can used as variant input, too. In this mode, --sample will not be used.
                        If multiple vcf files are provided, use "," to join the file names. For example, "--mutations file1.vcf,file2.vcf". A pattern match is also supported for input vcf, but quote is required to get it work. For example '--mutations "file*.vcf" ' 
                        
                        ''', default = '', required=False)
    parser.add_argument('-p','--protein', help = 'protein sequences in fasta format. It can be a gzip file. Only proteins in this file will be checked', default='')
    parser.add_argument('-t', '--threads', help='number of threads/CPUs to run the program. default, use 20 or all CPUs available, whichever is smaller', type=int, default=min(20, os.cpu_count()))
    parser.add_argument('-o', '--out', help='''output prefix, folder path could be included. Three or five files will be saved depending on the variant file format. Outputs include the annotation for mutated transcripts, the mutated or all protein sequences, two variant files from vcf. {out}.pergeno.aa_mutations.csv, {out}.pergeno.protein_all.fa, {out}.protein_changed.fa, {out}.vcf2mutation_1/2.tsv. default "perGeno" ''', default="perGeno")
    parser.add_argument('-a', '--datatype', help='''input datatype, could be GENCODE_GTF, GENCODE_GFF3, RefSeq, Ensembl_GTF or gtf. default "gtf". Ensembl_GFF3 is not supported. ''', default='gtf', type=str, choices=['GENCODE_GTF', 'GENCODE_GFF3','RefSeq','Ensembl_GTF','gtf'])
    parser.add_argument('-k','--protein_keyword', help='''field name in attribute column of gtf file to determine ids for proteins. default "auto", determine the protein_keyword based on datatype. "transcript_id" for GENCODE_GTF, "protein_id" for "RefSeq" and "Parent" for gtf and GENCODE_GFF3 ''', default='auto')
    parser.add_argument('-F', '--no_filter', help='default only keep variant with value "PASS" FILTER column of vcf file. if set, do not filter', action='store_true')
    parser.add_argument('-s', '--sample', help='''
                        sample name in the vcf/tsv to extract the variant information. default: None, extract the first sample in vcf file, or use all variants in the tsv file. 
                        For multiple samples, use "," to join the sample names. For example, "--sample sample1,sample2,sample3". 
                        To use all samples, use "--sample ALL_SAMPLES". 
                        To use all variants regardless where the variants from, use "--sample ALL_VARIANTS".
                        ''', default=None)
    parser.add_argument('-A','--all_chromosomes', help='default keep variant in chromosomes and ignore those in short fragments of the genome. if set, use all chromosomes including fragments when parsing the vcf file', action='store_true')
    parser.add_argument('-D','--download', help='''download could be 'GENCODE','RefSeq','Ensembl','Uniprot', 'CHM13'. If set, PrecisonProDB will try to download genome, gtf and protein files from the Internet. Download will be skipped if "--genome, --gtf, --protein, (--uniprot)" were all set. Settings from "--genome, --gtf, --protein, (--uniprot), --datatype" will not be used if the files were downloaded by PrecisonProDB. default "". Note, if --sqlite is set, will not download any files ''', default='', type=str, choices=['GENCODE','RefSeq','Ensembl','Uniprot','CHM13',''])
    parser.add_argument('-U','--uniprot', help='''uniprot protein sequences. If more than one file, use "," to join the files. default "". For example, "UP000005640_9606.fasta.gz", or "UP000005640_9606.fasta.gz,UP000005640_9606_additional.fasta" ''', default='', type=str)
    parser.add_argument('--uniprot_min_len', help='''minimum length required when matching uniprot sequences to proteins annotated in the genome. default 20 ''', default=20, type=int)
    parser.add_argument('--PEFF', help='If set, PEFF format file(s) will be generated. Default: do not generate PEFF file(s).', action='store_true')
    parser.add_argument('--keep_all', help='If set, do not delete files generated during the run', action='store_true')

    parser.add_argument('-S','--sqlite', help='''A path of sqlite file for re-use of annotation info. default '', do not use sqlite. The program will create a sqlite file if the file does not exist. If the file already exists, the program will use data stored in the file. It will cause error if the content in the sqlite file is not as expected. ''', default='', type=str)
    parser.add_argument('--info_field', help='fields to use in the INFO column of the vcf file to filter variants. Default None', default = None)
    parser.add_argument('--info_field_thres', help='threhold for the info field. Default None, do not filter any variants. If set "--info_filed AF --info_field_thres 0.01", only keep variants with AF >= 0.01', default = None)

    
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
    download = f.download
    files_uniprot = f.uniprot
    uniprot_min_len=f.uniprot_min_len
    keep_all = f.keep_all
    file_sqlite = f.sqlite
    print(f)

    
    time0 = time.time()

    # create workfolder if not exist
    workfolder = os.path.dirname(outprefix)
    if workfolder == '':
        workfolder = '.'
    if workfolder != '.':
        if not os.path.exists(workfolder):
            os.makedirs(workfolder)

    # download required files if download is set. will not download if file_sqlite is set. 
    download = download.upper()
    if download != '':
        print('-D --download is set to be', download, '\n')
        if download != 'UNIPROT':
            if file_genome != '' and file_gtf != '' and file_protein != '':
                print('download already finished for', download)
                pass
            else:
                files_downloaded = downloadHuman.download(download, workfolder=workfolder)
                file_genome, file_gtf, file_protein = files_downloaded
                if download == 'GENCODE':
                    datatype = 'GENCODE_GTF'
                elif download == 'REFSEQ':
                    datatype = 'RefSeq'
                elif download == 'ENSEMBL':
                    datatype = 'Ensembl_GTF'
                elif download == 'CHM13':
                    datatype = 'RefSeq'
                else:
                    print('download is', download, 'which is not possible...')
        else:
            if file_genome != '' and file_gtf != '' and file_protein != '' and files_uniprot != '':
                print('download already finished for', download)
                pass
            else:
                files_downloaded = downloadHuman.download(download, workfolder=workfolder)
                file_genome, file_gtf, file_protein, file_uniprot, file_uniprot_additional = files_downloaded
                files_uniprot = ','.join([file_uniprot, file_uniprot_additional])
                datatype = 'Ensembl_GTF'

    pattern = re.compile(r'(chr)?(\d+)-(\d+)-([A-Za-z]+)-([A-Za-z]+)')
    match = pattern.match(file_mutations)
    if match and (not os.path.exists(file_mutations)):
        if file_sqlite == '':
            print(f'file_mutations is a string {file_mutations} while file_sqlite is not provided. exit...')
            sys.exit()
            
    if individual == 'ALL_SAMPLES' or individual == "ALL_VARIANTS" or ',' in str(individual):
        if file_sqlite == '':
            print(f'sample is set to {individual}. In this case, --sqlite must be set. exit...')
            sys.exit()
    if ',' in file_mutations or '*' in file_mutations:
        if file_sqlite == '':
            print(f'mutations is set to {file_mutations}. In this case, --sqlite must be set. exit...')
            sys.exit()
    
    if file_sqlite == '':
        if file_mutations == '':
            print('file_sqlite not provided. no input mutation file is provided. exit...')
            sys.exit()
        if file_mutations.lower().endswith('.vcf') or file_mutations.lower().endswith('.vcf.gz'):
            print('variant file is a vcf file')
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
        else:
            print('variant file is a tsv file')
            pergeno = PerGeno(
                file_genome = file_genome, 
                file_gtf=file_gtf, 
                file_mutations = file_mutations, 
                file_protein=file_protein, 
                threads=threads, 
                outprefix=outprefix, 
                datatype=datatype, 
                protein_keyword=protein_keyword, 
                keep_all=keep_all
                )
            #print(pergeno.__dict__)
            pergeno.splitInputByChromosomes()
            #print(pergeno.__dict__)
            pergeno.runPerChom()
    else:
        # use Sqlite
        import PrecisionProDB_Sqlite
        print('using sqlite database to speed up')
        PrecisionProDB_Sqlite.main_PrecsionProDB_Sqlite(file_genome, file_gtf, file_mutations, file_protein, threads, outprefix, datatype, protein_keyword, filter_PASS, individual, chromosome_only, keep_all, file_sqlite, info_field = f.info_field, info_field_thres = f.info_field_thres)

    pattern = re.compile(r'(chr)?(\d+)-(\d+)-([A-Za-z]+)-([A-Za-z]+)')
    match = pattern.match(file_mutations)
    if match and (not os.path.exists(file_mutations)):
        if download == 'UNIPROT':
            print(f'file_mutations is a string {file_mutations}. running with UniProt is not supported. will not extract UniProt sequences or generate PEFF file!')
    else:
        # deal with uniprot
        if download == 'UNIPROT':
            if file_sqlite != '' and file_protein == '':
                print('extract all protein sequences from sqlite file')
                file_protein = outprefix + '.file_proteins_input_from_sqlite.fasta'
                buildSqlite.get_proteins_from_sqlite(file_sqlite, file_output = file_protein)

            print('try to extract Uniprot proteins from Ensembl models')
            import extractMutatedUniprot
            extractMutatedUniprot.extractMutatedUniprot(files_uniprot=files_uniprot, files_ref=file_protein, files_alt=outprefix + '.pergeno.protein_all.fa', outprefix=outprefix, length_min = uniprot_min_len)

        # generate PEFF output file
        if f.PEFF:
            import generatePEFFoutput
            generatePEFFoutput.generatePEFFoutput(file_protein = file_protein, file_mutation = outprefix + '.pergeno.aa_mutations.csv', file_out = outprefix + '.pergeno.protein_PEFF.fa', TEST=False, file_sqlite = file_sqlite)

            if download == 'UNIPROT':
                generatePEFFoutput.generateUniprotPEFFout(file_PEFF = outprefix + '.pergeno.protein_PEFF.fa', files_uniprot_ref = files_uniprot, file_uniprot_changed = outprefix + '.uniprot_changed.tsv', file_uniprot_out = outprefix + '.uniprot_PEFF.fa')


    print('PrecisionProDB finished! Total seconds:', time.time() - time0)

if __name__ == '__main__':
    main()