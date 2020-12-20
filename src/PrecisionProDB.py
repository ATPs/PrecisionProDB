from PrecisionProDB_core import PerGeno
from PrecisionProDB_vcf import runPerGenoVCF
import os
import downloadHuman
import time

description = '''
PrecisionProDB, personal proteogenomic tools which outputs a new reference protein based on the variants data. 
VCF/tsv file as the variant input
If the variant file is in tsv format, at least four columns are required in the header: chr, pos, ref, alt. Try to convert the file to proper format if you have a bed or other file.
'''
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-g','--genome', help = 'the reference genome sequence in fasta format. It can be a gzip file', default='')
    parser.add_argument('-f', '--gtf', help='gtf file with CDS and exon annotations. It can be a gzip file', default='')
    parser.add_argument('-m', '--mutations', help='''a file stores the variants. If the file ends with ".vcf" or ".vcf.gz", treat as vcf input. Otherwise, treat as TSV input ''', required=True)
    parser.add_argument('-p','--protein', help = 'protein sequences in fasta format. It can be a gzip file. Only proteins in this file will be checked', default='')
    parser.add_argument('-t', '--threads', help='number of threads/CPUs to run the program. default, use all CPUs available', type=int, default=os.cpu_count())
    parser.add_argument('-o', '--out', help='''output prefix. Five files will be saved, including the annotation for mutated transcripts, the mutated or all protein sequences, two variant files from vcf. {out}.pergeno.aa_mutations.csv, {out}.pergeno.protein_all.fa, {out}.protein_changed.fa, {out}.vcf2mutation_1/2.tsv. default "perGeno" ''', default="perGeno")
    parser.add_argument('-a', '--datatype', help='''input datatype, could be GENCODE_GTF, GENCODE_GFF3, RefSeq, Ensembl_GTF or gtf. default "gtf". Ensembl_GFF3 is not supported. ''', default='gtf', type=str, choices=['GENCODE_GTF', 'GENCODE_GFF3','RefSeq','Ensembl_GTF','gtf'])
    parser.add_argument('-k','--protein_keyword', help='''field name in attribute column of gtf file to determine ids for proteins. default "auto", determine the protein_keyword based on datatype. "transcript_id" for GENCODE_GTF, "protein_id" for "RefSeq" and "Parent" for gtf and GENCODE_GFF3 ''', default='auto')
    parser.add_argument('-F', '--no_filter', help='default only keep variant with value "PASS" FILTER column of vcf file. if set, do not filter', action='store_true')
    parser.add_argument('-s', '--sample', help='sample name in the vcf to extract the variant information. default: None, extract the first sample', default=None)
    parser.add_argument('-A','--all_chromosomes', help='default keep variant in chromosomes and ignore those in short fragments of the genome. if set, use all chromosomes including fragments when parsing the vcf file', action='store_true')
    parser.add_argument('-D','--download', help='''download could be 'GENCODE','RefSeq','Ensembl','Uniprot'. If set, perGeno will try to download genome, gtf and protein files from the Internet. Download will be skipped if "--genome, --gtf, --protein, (--uniprot)" were all set. Settings from "--genome, --gtf, --protein, (--uniprot), --datatype" will not be used if the files were downloaded by perGeno. default "".''', default='', type=str, choices=['GENCODE','RefSeq','Ensembl','Uniprot',''])
    parser.add_argument('-U','--uniprot', help='''uniprot protein sequences. If more than one file, use "," to join the files. default "". For example, "UP000005640_9606.fasta.gz", or "UP000005640_9606.fasta.gz,UP000005640_9606_additional.fasta" ''', default='', type=str)
    parser.add_argument('--uniprot_min_len', help='''minimum length required when matching uniprot sequences to proteins annotated in the genome. default 20 ''', default=20, type=int)
    
    
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
    print(f)
    
    time0 = time.time()

    # create workfolder if not exist
    workfolder = os.path.dirname(outprefix)
    if workfolder == '':
        workfolder = '.'
    if workfolder != '.':
        if not os.path.exists(workfolder):
            os.makedirs(workfolder)

    # download required files if download is set
    if download != '':
        print('-D --download is set to be', download, '\n')
        download = download.upper()
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



    if file_mutations.lower().endswith('.vcf') or file_mutations.lower().endswith('.vcf.gz'):
        print('variant file is a vcf file')
        runPerGenoVCF(file_genome = file_genome, file_gtf=file_gtf, file_mutations = file_mutations, file_protein=file_protein, threads=threads, outprefix=outprefix, datatype=datatype, protein_keyword=protein_keyword, filter_PASS=filter_PASS, individual=individual, chromosome_only=chromosome_only)
    else:
        print('variant file is a tsv file')
        pergeno = PerGeno(file_genome = file_genome, file_gtf=file_gtf, file_mutations = file_mutations, file_protein=file_protein, threads=threads, outprefix=outprefix, datatype=datatype, protein_keyword=protein_keyword)
        #print(pergeno.__dict__)
        pergeno.splitInputByChromosomes()
        #print(pergeno.__dict__)
        pergeno.runPerChom()
    
    # deal with uniprot
    if download == 'UNIPROT':
        print('try to extract Uniprot proteins from Ensembl models')
        import extractMutatedUniprot
        extractMutatedUniprot.extractMutatedUniprot(files_uniprot=files_uniprot, files_ref=file_protein, files_alt=outprefix + '.pergeno.protein_all.fa', outprefix=outprefix, length_min = uniprot_min_len)

    
    print('PrecisionProDB finished!, total seconds:', time.time() - time0)
