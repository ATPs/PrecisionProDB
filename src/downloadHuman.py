from ftplib import FTP
import os
import re

def getGENCODElatest():
    '''return urls of files from the latest release version of GENCODE
    '''
    ftp_url = 'ftp.ebi.ac.uk'
    folder = '/pub/databases/gencode/Gencode_human'
    ftp = FTP(ftp_url)
    ftp.login()
    folders = ftp.nlst(folder)
    folders = [e for e in folders if '/release_' in e]
    versions = [e.split('/release_')[1] for e in folders]
    v2 = []
    for e in versions:
        try:
            i = [int(e), '']
        except:
            i = [int(e[:-1]), e[-1]]
        v2.append(i)
    
    v2 = sorted(v2)
    latest = v2[-1]
    latest = str(latest[0]) + latest[1]
    print('GENCODE, latest release is', latest)
    
    folder_latest = [e for e in folders if '/release_'+latest in e][0]
    files_latest = ftp.nlst(folder_latest)

    url_genome = 'ftp://' + ftp_url + [e for e in files_latest if re.findall('GRCh38.p\\d*.genome.fa.gz',e)][0]
    url_GTF = 'ftp://' + ftp_url + [e for e in files_latest if re.findall('gencode.v\\d*.chr_patch_hapl_scaff.annotation.gtf.gz',e)][0]
    url_protein = 'ftp://' + ftp_url + [e for e in files_latest if re.findall('gencode.v\\d*.pc_translations.fa.gz',e)][0]

    return url_genome, url_GTF, url_protein


def getRefSeqLatest():
    '''return urls of files from the latest release version of RefSeq
    '''
    ftp_url = 'ftp.ncbi.nlm.nih.gov'
    folder = '/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/current'
    ftp = FTP(ftp_url)
    ftp.login()
    folders = ftp.nlst(folder)

    # get the only folder in folders
    for file in folders:
        try:
            ftp.size(file)
        except:
            break
    
    folder_release = file
    version = os.path.basename(folder_release)
    print('RefSeq, the latest release is', version)

    files_latest = ftp.nlst(folder_release)

    url_genome = 'ftp://' + ftp_url + [e for e in files_latest if re.findall(f'/{version}_genomic.fna.gz',e)][0]
    url_GTF = 'ftp://' + ftp_url + [e for e in files_latest if re.findall(f'/{version}_genomic.gtf.gz',e)][0]
    url_protein = 'ftp://' + ftp_url + [e for e in files_latest if re.findall(f'/{version}_protein.faa.gz',e)][0]

    return url_genome, url_GTF, url_protein

def getEnsemblLatest():
    '''return urls of files from the latest release version of Ensembl
    '''
    ftp_url = 'ftp.ensembl.org'
    ftp = FTP(ftp_url)
    ftp.login()
    folder_genome = '/pub/current_fasta/homo_sapiens/dna/'
    folder_GTF = '/pub/current_gtf/homo_sapiens/'
    folder_protein = '/pub/current_fasta/homo_sapiens/pep/'

    files_genome = ftp.nlst(folder_genome)
    files_GTF = ftp.nlst(folder_GTF)
    files_protein = ftp.nlst(folder_protein)

    url_genome = 'ftp://' + ftp_url + [e for e in files_genome if re.findall(f'/Homo_sapiens.GRCh\\d*.dna.primary_assembly.fa.gz',e)][0]
    url_GTF = 'ftp://' + ftp_url + [e for e in files_GTF if re.findall(f'Homo_sapiens.GRCh\\d*.\\d*.gtf.gz',e)][0]
    url_protein = 'ftp://' + ftp_url + [e for e in files_protein if re.findall(f'Homo_sapiens.GRCh\\d*.pep.all.fa.gz',e)][0]

    version = os.path.basename(url_GTF).split('.')[2]
    print('Ensembl, the latest release is', version)

    return url_genome, url_GTF, url_protein

def downloadFile2folder(url, workfolder):
    '''
    download the file from url to workfolder. print url if cannot download
    '''
    try:
        # basename = os.path.basename(url)
        # if os.path.exists(os.path.join(workfolder, basename)):
        #     print(basename, 'already downloaded in folder', workfolder, 'and will skip')
        #     return 0
        os.system(f'cd {workfolder} && wget -c --retry-connrefused --waitretry=2 --read-timeout=40 --timeout=15 -t 20 {url}')
    except:
        print(url, 'cannot be downloaded, try to download yourself!')

def download(datatype, workfolder='.'):
    '''download genome, GTF, protein sequence of {datatype} and store files in {workfolder}
    datatype can be RefSeq, GENCODE or Ensembl
    default workfolder is current folder
    '''
    datatype = datatype.upper()
    if datatype == 'GENCODE':
        url_genome, url_GTF, url_protein = getGENCODElatest()
    elif datatype == 'REFSEQ':
        url_genome, url_GTF, url_protein = getRefSeqLatest()
    elif datatype == 'ENSEMBL':
        url_genome, url_GTF, url_protein = getEnsemblLatest()
    elif datatype == 'UNIPROT':
        url_genome, url_GTF, url_protein = getEnsemblLatest()
    else:
        print('other datatype is not supported now.')
        exit(0)
    
    if not os.path.exists(workfolder):
        os.makedirs(workfolder)

    downloadFile2folder(url_genome, workfolder)
    downloadFile2folder(url_GTF, workfolder)
    downloadFile2folder(url_protein, workfolder)
    

    if datatype == 'UNIPROT':
        url_uniprot = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640_9606.fasta.gz'
        downloadFile2folder(url_uniprot, workfolder)
        url_uniprot_additional = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640_9606_additional.fasta.gz'
        downloadFile2folder(url_uniprot_additional, workfolder)



    print('finished! Files saved in', workfolder)
    # return path of downloaded files
    if datatype != 'UNIPROT':
        urls = [url_genome, url_GTF, url_protein]
    else:
        urls = [url_genome, url_GTF, url_protein, url_uniprot, url_uniprot_additional]
    files_download = [os.path.join(workfolder, os.path.basename(e)) for e in urls]
    return files_download


description = '''
download the latest human gene models from RefSeq, GENCODE or Ensembl to run perGeno
If datatype is "Uniprot", Ensembl and Uniprot human sequences (UP000005640_9606, UP000005640_9606_additional) will be downloaded
'''

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-d','--datatype', help = 'RefSeq, GENCODE, Ensembl or Uniprot to download.', required=True)
    parser.add_argument('-o', '--out', help='output folder for the downloaded files. default: ".",current folder', default='.')
    f = parser.parse_args()
    
    datatype = f.datatype
    workfolder = f.out
    download(datatype=datatype, workfolder=workfolder)

