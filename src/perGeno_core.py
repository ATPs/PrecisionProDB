import pandas as pd
from Bio import SeqIO
import gzip
import os
import time


def openFile(filename):
    '''open text or gzipped text file
    '''
    if filename.endswith('.gz'):
        return gzip.open(filename,'rt')
    return open(filename)


def iterTxtWithComment(filename, comment = '#'):
    '''filename is a text file or a gzipped text file
    iterate lines, excluding empty lines and lines startswith comment.
    newline symbol removed.
    default comment string is '#'
    '''
    fo = openFile(filename)
    for line in fo:
        if line.startswith('#'):
            continue
        line = line.strip()
        if len(line) == 0:
            continue
        yield line


class PerGeno(object):
    """
    Define an instance of PerGeno analysis
    contains all data, all parameters and provides the main methods for permoing analysis
    """
    def __init__(
                    self,
                    file_genome,
                    file_gtf,
                    file_mutations,
                    file_protein,
                    threads = os.cpu_count(),
                    outprefix = 'perGeno',
                    datatype = 'gtf',
                    protein_keyword = 'auto'
                ):
        # variable from input
        self.file_genome = file_genome # genome file location
        self.file_gtf = file_gtf # gtf file location
        self.file_mutations = file_mutations # mutation file location
        self.file_protein = file_protein # protein file location
        self.threads = threads # threads to use when running program
        self.outprefix = outprefix # where to store the results
        self.datatype = datatype # input datatype, could be GENCODE_GTF, GENCODE_GFF3, RefSeq or gtf
        self.protein_keyword = protein_keyword #protein_keyword, field name in attribute column of gtf file to determine ids for proteins. default "auto", determine the protein_keyword based on datatype. "transcript_id" for GENCODE_GTF, "protein_id" for "RefSeq" and "Parent" for gtf and GENCODE_GFF3

        # determine datatype
        self.datatype = self.determineInputType()
        
        # protein_keyword, field name in attribute column of gtf file
        if protein_keyword != 'auto':
            self.protein_keyword = protein_keyword
        else:
            if self.datatype == 'GENCODE_GTF':
                self.protein_keyword = 'transcript_id'
            elif self.datatype == 'GENCODE_GFF3':
                self.protein_keyword = 'Parent'
            elif self.datatype == 'RefSeq':
                self.protein_keyword = 'protein_id'
            else:
                self.protein_keyword = 'Parent'
        
        # folder to store temporary data
        self.tempfolder = self.outprefix +'_temp'
        
        # prepare output folder
        if not os.path.exists(self.tempfolder):
            os.makedirs(self.tempfolder)
        
        # some other chromosomes info to store
        self.chromosomes_genome = None
        self.chromosomes_genome_description = None
        self.chromosomes_mutation = None
        self.chromosomes_gtf = None

    def determineInputType(self):
        '''determine datatype based on gtf/gff file
        '''
        datatype = self.datatype
        file_gtf = self.file_gtf

        # datatype is already provided
        if datatype.upper() == 'GENCODE_GTF':
            datatype = 'GENCODE_GTF'
            print('input data is provided as from GENCODE_GTF')
            return datatype
        if datatype.upper() == 'GENCODE_GFF3':
            datatype = 'GENCODE_GFF3'
            print('input data is provided as from GENCODE_GFF3')
            return datatype
        if datatype.upper() == 'REFSEQ':
            datatype = 'RefSeq'
            print('input data is provided as from RefSeq')
            return datatype

        f = openFile(file_gtf)
        annotations = ''
        for line in f:
            if line.startswith('#'):
                annotations += line
            else:
                break
        
        annotations = annotations.lower()
        # check if it is gencode
        if 'gencode' in annotations:
            if 'format: gff3' in annotations:
                datatype = 'GENCODE_GFF3'
                print('input data is inferred as from GENCODE_GFF3')
                return datatype
            elif 'format: gtf' in annotations:
                datatype = 'GENCODE_GTF'
                print('input data is inferred as from GENCODE_GTF')
                return datatype
            else:
                datatype = 'GENCODE_GTF'
                print('input data is inferred as from GENCODE, but cannot determine from GFF3 or GTF. treat as GENCODE_GTF')
                return datatype
        if 'NCBI Homo sapiens'.lower() in annotations:
            datatype = 'RefSeq'
            print('input data is from NCBI Homo sapiens RefSeq')
            return datatype
        
        datatype = 'gtf'
        print('input data is general gtf')
        return datatype


    def getProteinIDsForEachChromosome(self):
        '''split proteins based on chromosomes.
        GENCODE, use transcript_id as protein name
        for RefSeq or gtf, use protein_id as protein name
        '''
        datatype = self.datatype
        file_proteins = self.file_protein
        protein_keyword = self.protein_keyword
        file_gtf = self.file_gtf

        protein_ids = set([e.id for e in SeqIO.parse(openFile(file_proteins), 'fasta')])
        keyword = protein_keyword
        dc_protein2chr = {}
        if 'GENCODE' in datatype:
            protein_ids = set([e.split('|')[1] for e in protein_ids])
        
        for line in iterTxtWithComment(file_gtf):
            es = line.split('\t')
            chromosome, anno = es[0], es[-1]
            elements = anno.strip(';').split(';')
            elements = [e.strip() for e in elements]
            for element in elements:
                if keyword in element:
                    if '=' in element:
                        protein_id = element.split('=')[1]
                    else:
                        protein_id = element.split(' ', maxsplit=1)[1].strip('"')
                    if protein_id in protein_ids:
                        dc_protein2chr[protein_id] = chromosome
        
        if len(dc_protein2chr) == 0:# if "Parent not in attribute field of gtf"
            keyword = 'protein_id'
        
        for line in iterTxtWithComment(file_gtf):
            es = line.split('\t')
            chromosome, anno = es[0], es[-1]
            elements = anno.strip(';').split(';')
            elements = [e.strip() for e in elements]
            for element in elements:
                if keyword in element:
                    if '=' in element:
                        protein_id = element.split('=')[1]
                    else:
                        protein_id = element.split(' ', maxsplit=1)[1].strip('"')
                    if protein_id in protein_ids:
                        dc_protein2chr[protein_id] = chromosome
        
        return dc_protein2chr

    def splitGenomeByChromosomes(self):
        '''split the genome based on chromosomes
        '''
        file_genome = self.file_genome
        tempfolder = self.tempfolder

        chromosomes_genome = []
        chromosomes_genome_description = []
        for e in SeqIO.parse(openFile(file_genome),'fasta'):
            chromosomes_genome.append(e.id)
            chromosomes_genome_description.append(e.description)
            tf = os.path.join(tempfolder,e.id + '.genome.fa')
            open(tf, 'w').write('>{}\n{}\n'.format(e.id, str(e.seq).strip('*')))
        print('finish splitting the genome file')
        return chromosomes_genome, chromosomes_genome_description

    def splitProteinByChromosomes(self,dc_protein2chr):
        '''split the protein file based on chromosomes
        '''
        datatype = self.datatype
        tempfolder = self.tempfolder
        file_protein = self.file_protein

        chromosomes_protein = []
        dc_files = {}
        tf_proteins_not_in_genome = os.path.join(tempfolder, 'not_in_genome.proteins.fa')
        tf_proteins_not_in_genome = open(tf_proteins_not_in_genome,'w')
        for seq in SeqIO.parse(openFile(file_protein),'fasta'):
            if 'GENCODE' in datatype:# use transcript_id as protein_id
                protein_id = seq.id.split('|')[1]
            else:
                protein_id = seq.id
            if protein_id in dc_protein2chr:
                chromosome = dc_protein2chr[protein_id]
                if chromosome not in dc_files:
                    tf = os.path.join(tempfolder, chromosome + '.proteins.fa')
                    dc_files[chromosome] = open(tf,'w')
                dc_files[chromosome].write('>{}\t{}\n{}\n'.format(protein_id, seq.description, str(seq.seq).strip('*')))# remove "*" which represent stop codons
            else:
                tf_proteins_not_in_genome.write('>{}\t{}\n{}\n'.format(protein_id, seq.description, str(seq.seq).strip('*')))# remove "*" which represent stop codons
        for chromosome in dc_files:
            chromosomes_protein.append(chromosome)
            dc_files[chromosome].close()
        print('finish splitting the protein file')
        tf_proteins_not_in_genome.close()
        return chromosomes_protein

    def splitMutationByChromosome(self):
        '''split mutation file based on chromosomes
        '''
        tempfolder = self.tempfolder
        file_mutations = self.file_mutations
        chromosomes_genome = self.chromosomes_genome
        chromosomes_genome_description = self.chromosomes_genome_description

        df_mutations = pd.read_csv(file_mutations, sep='\t',low_memory=False)
        df_mutations['chr'] = df_mutations['chr'].astype(str)
        chromosomes_mutation = []

        for k,v in df_mutations.groupby('chr'):
            v = v.copy()
            chromosomes_mutation.append(k)
            if k in chromosomes_genome:
                tf = os.path.join(tempfolder, k + '.mutation.tsv')
            elif 'chr' + k in chromosomes_genome:
                print('add "chr" to chromosome ' + k +' in mutation file')
                tf = os.path.join(tempfolder, 'chr' + k + '.mutation.tsv')
                v['chr'] = v['chr'].apply(lambda x:'chr' + str(x))
            else:
                print('chromosomes in mutation file is different from the genome. try to solve that. This is usually True if datatype is RefSeq')
                if k.startswith('chr'):
                    k = k[3:]
                if k == 'M':
                    for e in chromosomes_genome_description:
                        e1, e2 = e.split(' ', maxsplit=1)
                        if 'Homo sapiens mitochondrion, complete genome' in e2:
                            k_new = e1
                            print(f'    mutation chromosome change {k} to {k_new}')
                            break
                else:
                    for e in chromosomes_genome_description:
                        e1, e2 = e.split(' ', maxsplit=1)
                        if f'Homo sapiens chromosome {k},' in e2:
                            k_new = e1
                            print(f'    mutation chromosome change {k} to {k_new}')
                            break
                    else:
                        print('chromosome in mutation file', k, 'cannot find corresponding chromosome in genome file. please check')
                        k_new = k
                tf = os.path.join(tempfolder, k_new + '.mutation.tsv')

            v.to_csv(tf, sep='\t',index=None)
        print('finish splitting the mutation file')
        return chromosomes_mutation

    def splitGtfByChromosomes(self,dc_protein2chr):
        '''split gtf file based on chromosome. only keep proteins in file_protein
        '''
        file_gtf = self.file_gtf
        tempfolder = self.tempfolder
        protein_keyword = self.protein_keyword
        
        chromosomes_gtf = []
        dc_files = {}
        features_keep = ['CDS', 'exon']
        
        for line in iterTxtWithComment(file_gtf):
            es = line.strip().split('\t')
            chromosome = es[0]
            feature = es[2]
            elements = es[-1].strip(';').split(';')
            elements = [e.strip() for e in elements]
            protein_ID = False
            for element in elements:
                if protein_keyword in element:
                    if '=' in element:
                        protein_ID = element.split('=')[1]
                    else:
                        protein_ID = element.split(' ', maxsplit=1)[1].strip('"')
            if chromosome not in dc_files:
                tf = os.path.join(tempfolder, chromosome + '.gtf')
                dc_files[chromosome] = open(tf,'w')
            
            if feature in features_keep and protein_ID:
                if protein_ID in dc_protein2chr:
                    es[-1] = protein_ID
                    dc_files[chromosome].write('\t'.join(es) + '\n')
        for chromosome in dc_files:
            chromosomes_gtf.append(chromosome)
            dc_files[chromosome].close()
        print('finish splitting the gtf file')
        return chromosomes_gtf


    def splitInputByChromosomes(self):
        '''split the genome, mutation, gtf and protein file based on chromosomes
        '''
        print('split the genome, mutation, gtf and protein file based on chromosomes')
        
        # split genome
        self.chromosomes_genome,self.chromosomes_genome_description = self.splitGenomeByChromosomes()

        # split protein sequences
        dc_protein2chr = self.getProteinIDsForEachChromosome()
        print('finish assign chromosome for each protein')
        self.chromosomes_protein = self.splitProteinByChromosomes(dc_protein2chr)

        # split mutation. if the chromosome column 'chr' is not the same as chromosomes_genome, try to fix that.
        self.chromosomes_mutation = self.splitMutationByChromosome()
        
        # split gtf
        self.chromosomes_gtf = self.splitGtfByChromosomes(dc_protein2chr)

        self.chromosomes = self.chromosomes_protein
        print('number of chromosomes with proteins:', len(self.chromosomes))
        return None
    


description = '''PerGeno, personal proteogenomic tools which outputs a new reference protein based on the variants data
'''

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-g','--genome', help = 'the reference genome sequence in fasta format. It can be a gzip file', required=True)
    parser.add_argument('-f', '--gtf', help='gtf file with CDS and exon annotations', required=True)
    parser.add_argument('-m', '--mutations', help='a file stores the variants', required=True)
    parser.add_argument('-p','--protein', help = 'protein sequences in fasta format. It can be a gzip file. Only proteins in this file will be checked', required=True)
    parser.add_argument('-t', '--threads', help='number of threads/CPUs to run the program. default, use all CPUs available', type=int, default=os.cpu_count())
    parser.add_argument('-o', '--out', help='''output prefix. two file will be output. One is the annotation for mutated transcripts, one is the protein sequences. {out}_mutations.csv, {out}_protein.fa. default "perGeno" ''', default="perGeno")
    parser.add_argument('-a', '--datatype', help='''input datatype, could be GENCODE, RefSeq or gtf. default "gtf" ''', default='gtf', type=str, choices=['GENCODE_GTF', 'GENCODE_GFF3','RefSeq','gtf'])
    parser.add_argument('-k','--protein_keyword', help='''field name in attribute column of gtf file to determine ids for proteins. default "auto", determine the protein_keyword based on datatype. "transcript_id" for GENCODE_GTF, "protein_id" for "RefSeq" and "Parent" for gtf and GENCODE_GFF3 ''', default='auto')
    f = parser.parse_args()
    
    file_genome = f.genome
    file_gtf = f.gtf
    file_mutations = f.mutations
    file_protein = f.protein
    threads = f.threads
    outprefix = f.out
    datatype = f.datatype
    protein_keyword = f.protein_keyword
    
    pergeno = PerGeno(file_genome = file_genome, file_gtf=file_gtf, file_mutations = file_mutations, file_protein=file_protein, threads=threads, outprefix=outprefix, datatype=datatype, protein_keyword=protein_keyword)
    print(pergeno.__dict__)
    pergeno.splitInputByChromosomes()
    print(pergeno.__dict__)
