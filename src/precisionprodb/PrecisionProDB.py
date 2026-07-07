import os
import time
import sys
import re

if __package__:
    from .args import add_argument_set
    from . import buildSqlite
    from . import downloadHuman
    from .PrecisionProDB_core import PerGeno
    from .PrecisionProDB_vcf import runPerGenoVCF
    from .vcf2mutation import is_manifest_file
else:
    from args import add_argument_set
    import buildSqlite
    import downloadHuman
    from PrecisionProDB_core import PerGeno
    from PrecisionProDB_vcf import runPerGenoVCF
    from vcf2mutation import is_manifest_file

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


def build_parser():
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-v', '--version', action='version', version=f'PrecisionProDB {get_version()}')
    add_argument_set(
        parser,
        'reference_inputs',
        'variant_input_detailed',
        'runtime_output',
        'annotation',
        'vcf_selection_mixed',
        'download_uniprot',
        'cleanup',
        'sqlite',
        'vcf_info_filters',
    )
    return parser


def run_from_args(f):
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
    uniprot_min_len = f.uniprot_min_len
    keep_all = f.keep_all
    file_sqlite = f.sqlite
    print(f)
    if file_sqlite == '':
        file_sqlite = outprefix + '.sqlite'
    elif file_sqlite == 'NONE':
        file_sqlite = ''

    
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
            reference_files_ready = file_genome != '' and file_gtf != '' and file_protein != ''
            sqlite_reference_ready = file_sqlite != '' and os.path.exists(file_sqlite)
            if files_uniprot != '' and (reference_files_ready or sqlite_reference_ready):
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
    if ',' in file_mutations or '*' in file_mutations or is_manifest_file(file_mutations):
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
        if __package__:
            from . import PrecisionProDB_Sqlite
        else:
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
            if __package__:
                from . import extractMutatedUniprot
            else:
                import extractMutatedUniprot
            extractMutatedUniprot.extractMutatedUniprot(files_uniprot=files_uniprot, files_ref=file_protein, files_alt=outprefix + '.pergeno.protein_all.fa', outprefix=outprefix, length_min = uniprot_min_len)

        # generate PEFF output file
        if f.PEFF:
            if __package__:
                from . import generatePEFFoutput
            else:
                import generatePEFFoutput
            generatePEFFoutput.generatePEFFoutput(file_protein = file_protein, file_mutation = outprefix + '.pergeno.aa_mutations.csv', file_out = outprefix + '.pergeno.protein_PEFF.fa', TEST=False, file_sqlite = file_sqlite)

            if download == 'UNIPROT':
                generatePEFFoutput.generateUniprotPEFFout(file_PEFF = outprefix + '.pergeno.protein_PEFF.fa', files_uniprot_ref = files_uniprot, file_uniprot_changed = outprefix + '.uniprot_changed.tsv', file_uniprot_out = outprefix + '.uniprot_PEFF.fa')


    print('PrecisionProDB finished! Total seconds:', time.time() - time0)


def main(argv=None):
    parser = build_parser()
    run_from_args(parser.parse_args(argv))


if __name__ == '__main__':
    main()
