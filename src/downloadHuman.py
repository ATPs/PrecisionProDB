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
    ls_folders = []
    for file in folders:
        try:
            ftp.size(file)
        except:
            ls_folders.append(file)
    
    ls_url_genome = []
    ls_url_GTF = []
    ls_url_protein = []

    for folder_release in ls_folders:
        version = os.path.basename(folder_release)
        print('RefSeq, the latest release is', version)

        files_latest = ftp.nlst(folder_release)

        ls_url_genome += ['https://' + ftp_url + e for e in files_latest if re.findall(f'/.*_genomic.fna.gz',e)]
        ls_url_GTF += ['https://' + ftp_url + e for e in files_latest if re.findall(f'/.*_genomic.gtf.gz',e)]
        ls_url_protein += ['https://' + ftp_url + e for e in files_latest if re.findall(f'/.*_protein.faa.gz',e)]

    url_genome = [i for i in ls_url_genome if 'T2T-CHM13' not in i and 'cds_' not in i and 'rna_' not in i][0]
    url_GTF = [i for i in ls_url_GTF if 'T2T-CHM13' not in i][0]
    url_protein = [i for i in ls_url_protein if 'T2T-CHM13' not in i][0]
    return url_genome, url_GTF, url_protein


def getRefSeqLatestCHM13():
    '''return urls of files from the latest release version of RefSeq
    '''
    ftp_url = 'ftp.ncbi.nlm.nih.gov'
    folder = '/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/current'
    ftp = FTP(ftp_url)
    ftp.login()
    folders = ftp.nlst(folder)

    # get the only folder in folders
    ls_folders = []
    for file in folders:
        try:
            ftp.size(file)
        except:
            ls_folders.append(file)
    
    ls_url_genome = []
    ls_url_GTF = []
    ls_url_protein = []

    for folder_release in ls_folders:
        version = os.path.basename(folder_release)
        print('RefSeq, the latest release is', version)

        files_latest = ftp.nlst(folder_release)

        ls_url_genome += ['https://' + ftp_url + e for e in files_latest if re.findall(f'/.*_genomic.fna.gz',e)]
        ls_url_GTF += ['https://' + ftp_url + e for e in files_latest if re.findall(f'/.*_genomic.gtf.gz',e)]
        ls_url_protein += ['https://' + ftp_url + e for e in files_latest if re.findall(f'/.*_protein.faa.gz',e)]

    url_genome = [i for i in ls_url_genome if 'T2T-CHM13' in i and 'cds_' not in i and 'rna_' not in i][0]
    url_GTF = [i for i in ls_url_GTF if 'T2T-CHM13' in i][0]
    url_protein = [i for i in ls_url_protein if 'T2T-CHM13' in i][0]
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

def getUniprotLatest():
    '''return urls of files of latest Uniprot human sequences
    '''
    ftp_url = 'ftp.uniprot.org'
    ftp = FTP(ftp_url)
    ftp.login()
    folder = '/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/'
    files = ftp.nlst(folder)

    files = [e for e in files if 'UP000005640' in e]
    files2 = []
    for file in files:
        # get the only folder in folders
        try:
            ftp.size(file)
            files2.append(file)
        except:
            print(file, 'is a folder')
            for f2 in ftp.nlst(file):
                try:
                    ftp.size(f2)
                    files2.append(f2)
                except:
                    print(f2, 'is a folder')

    url_uniprot = [e for e in files2 if 'UP000005640_9606.fasta' in e][0]
    url_uniprot_additional = [e for e in files2 if 'UP000005640_9606_additional.fasta' in e][0]
    url_uniprot = 'ftp://' + ftp_url + url_uniprot
    url_uniprot_additional = 'ftp://' + ftp_url + url_uniprot_additional

    return url_uniprot, url_uniprot_additional


def downloadFile2folder(url, workfolder):
    '''
    download the file from url to workfolder. print url if cannot download
    '''
    try:
        # basename = os.path.basename(url)
        # if os.path.exists(os.path.join(workfolder, basename)):
        #     print(basename, 'already downloaded in folder', workfolder, 'and will skip')
        #     return 0
        a = os.system(f'cd {workfolder} && wget -c --retry-connrefused --waitretry=2 --read-timeout=40 --timeout=15 -t 20 {url}')
        # if not linux system, try to download with python lib
        if a != 0:
            import urllib 
            urllib.request.urlretrieve(url, os.path.join(workfolder, os.path.basename(url)))
    except:
        print(url, 'cannot be downloaded, try to download yourself! Check the https://github.com/ATPs/PrecisionProDB for information. Use the latest version or contact the authors!')
        exit()
    # check if the file is downloaded sucessfully
    basename = os.path.basename(url)
    if not os.path.exists(os.path.join(workfolder, basename)):
        print(url, 'cannot be downloaded, try to download yourself! Check the https://github.com/ATPs/PrecisionProDB for information. Use the latest version or contact the authors!')
        exit()

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
    elif datatype == 'CHM13':
        url_genome, url_GTF, url_protein = getRefSeqLatestCHM13()
    else:
        print('other datatype is not supported now.')
        exit(0)
    
    if not os.path.exists(workfolder):
        os.makedirs(workfolder)

    downloadFile2folder(url_genome, workfolder)
    downloadFile2folder(url_GTF, workfolder)
    downloadFile2folder(url_protein, workfolder)
    

    if datatype == 'UNIPROT':
        url_uniprot, url_uniprot_additional = getUniprotLatest()
        downloadFile2folder(url_uniprot, workfolder)
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
download the latest human gene models from RefSeq, GENCODE, Ensembl or UniProt to run PrecisionProDB.
If datatype is "Uniprot", Ensembl and UniProt human sequences (UP000005640_9606, UP000005640_9606_additional) will be downloaded.
If datatype is 'CHM13', the RefSeq CHM13 file will be downloaded
'''
def main():
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-d','--datatype', help = 'RefSeq, CHM13, GENCODE, Ensembl or Uniprot to download.', required=True)
    parser.add_argument('-o', '--out', help='output folder for the downloaded files. default: ".",current folder', default='.')
    f = parser.parse_args()
    
    datatype = f.datatype
    workfolder = f.out
    download(datatype=datatype, workfolder=workfolder)

if __name__ == '__main__':
    main()