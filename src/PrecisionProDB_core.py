import pandas as pd
from Bio import SeqIO
import gzip
import os
import time
import pickle
from perChrom import PerChrom
import shutil
import re

def get_k_new(k, chromosomes_genome, chromosomes_genome_description):
    '''k is chromosome name. return chromosome name based on chromosomes_genome, chromosomes_genome_description
    deal with RefSeq IDs
    '''
    if k.startswith('chr'):
        k = k[3:]
    # Step 1: Check if k or 'chr' + k is in chromosomes_genome
    if k in chromosomes_genome:
        return k
    elif 'chr' + k in chromosomes_genome:
        return 'chr' + k

    # Step 2: Special handling for mitochondrial chromosome 'M'
    if k == 'M':
        for e in chromosomes_genome_description:
            e1, e2 = e.split(' ', maxsplit=1)
            if 'mitochondrion' in e2.lower():
                return e1  # Return corresponding genome name

    # Step 3: General search for other chromosomes in chromosomes_genome_description
    for e in chromosomes_genome_description:
        e1, e2 = e.split(' ', maxsplit=1)
        if f'chromosome {k},' in e2:
            return e1  # Return corresponding genome name

    # Step 4: Fallback - if no match, return k as-is
    return k


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

def parseGtfAttributes(anno):
    '''anno is the last field of a gtf file or gff3 file, return a dict by parsing the contents
    '''
    elements = anno.strip().strip(';').split(';')
    elements = [e.strip() for e in elements]
    # tdc to store values in anno
    tdc = {}
    for element in elements:
        if '=' in element:
            k,v = element.split('=')
        else:
            k,v = element.split(' ', maxsplit=1)
        v = v.strip().strip('"')
        tdc[k] = v
    
    return tdc


class PerGeno(object):
    """
    Define an instance of PerGeno analysis
    contains all data, all parameters and provides the main methods for performing analysis
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
                    protein_keyword = 'auto',
                    keep_all = False
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
        self.keep_all = keep_all

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
            elif self.datatype == 'RefSeq' or self.datatype == 'Ensembl_GTF':
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
        if datatype.upper() == 'ENSEMBL_GTF':
            datatype = 'Ensembl_GTF'
            print('input data is provided as from Ensembl_GTF. Note: for some Ensembl models, the protein sequences are not translated based on GTF annotation. Those sequences will be ignored.')
            return datatype
        
        if not os.path.exists(file_gtf):
            return datatype
        print('datatype is not provided. Try to infer datatype. Note: Ensembl datatype cannot be inferred.')

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
            print('input data is inferred as from NCBI Homo sapiens RefSeq')
            return datatype
        

        # try to check if Ensembl
        line1 = f.readline()
        if re.findall('''\tgene_id "ENSG\\d*";''', line1):
            datatype = 'Ensembl_GTF'
            print('input data is inferred as from Ensembl_GTF. Use with caution!')
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

        # check if already exist protein2chr
        file_protein2chr = os.path.join(self.tempfolder,'protein2chr.pickle')
        if os.path.exists(file_protein2chr):
            print(f'protein-chromosome is already determined in file {file_protein2chr}, use previous result')
            dc_protein2chr = pickle.load(open(file_protein2chr,'rb'))
            return dc_protein2chr

        protein_ids = set([e.id for e in SeqIO.parse(openFile(file_proteins), 'fasta')])
        keyword = protein_keyword
        dc_protein2chr = {}
        dc_protein2multipleChr = {}
        if 'GENCODE' in datatype:
            protein_ids = set([e.split('|')[1] for e in protein_ids])
        if datatype == 'Ensembl_GTF':
            protein_ids = set([e.split('.')[0] for e in protein_ids])
        
        for line in iterTxtWithComment(file_gtf):
            es = line.split('\t')
            chromosome, feature, anno = es[0], es[2], es[-1]
            if feature != 'CDS':
                continue
            
            tdc = parseGtfAttributes(anno)
            if keyword in tdc:
                protein_id = tdc[keyword]
                if protein_id in protein_ids:
                    if protein_id not in dc_protein2chr:
                        dc_protein2chr[protein_id] = chromosome
                    else:
                        if chromosome != dc_protein2chr[protein_id]:
                            if protein_id not in dc_protein2multipleChr:
                                dc_protein2multipleChr[protein_id] = {dc_protein2chr[protein_id]}
                            dc_protein2multipleChr[protein_id].add(chromosome)
        
        if len(dc_protein2multipleChr) != 0:
            print(len(dc_protein2multipleChr), 'proteins exist in more than one chromosome, use the first chromosome in the gtf file. Choose a different protein_keyword if this is not the desired behavior. Current protein_keyword is', protein_keyword, '\nThose proteins in multiple chromosomes are:')
            for k,v in dc_protein2multipleChr.items():
                print(k,v)
        
        if len(dc_protein2chr) == 0:
            print('cannot determine chromosome for proteins based on protein_keyword:', protein_keyword)

        # write dc_protein2chr, as a checkpoint file
        with open(file_protein2chr,'wb') as f:
            pickle.dump(dc_protein2chr, f, protocol=4)

        return dc_protein2chr

    def splitGenomeByChromosomes(self):
        '''split the genome based on chromosomes
        '''
        file_genome = self.file_genome
        tempfolder = self.tempfolder

        # check if genome is already splitted
        file_splitGenomeFinished = os.path.join(self.tempfolder,'splitGenomeFinished.pickle')
        if os.path.exists(file_splitGenomeFinished):
            print(f'chromosome is already splitted. info stored in {file_splitGenomeFinished}, use previous result')
            chromosomes_genome, chromosomes_genome_description = pickle.load(open(file_splitGenomeFinished,'rb'))
            return chromosomes_genome, chromosomes_genome_description

        chromosomes_genome = []
        chromosomes_genome_description = []
        for e in SeqIO.parse(openFile(file_genome),'fasta'):
            chromosomes_genome.append(e.id)
            chromosomes_genome_description.append(e.description)
            tf = os.path.join(tempfolder,e.id + '.genome.fa')
            open(tf, 'w').write('>{}\n{}\n'.format(e.description, str(e.seq).strip('*')))
        print('finish splitting the genome file')

        # add checkpoint file
        with open(file_splitGenomeFinished,'wb') as f:
            pickle.dump([chromosomes_genome, chromosomes_genome_description], f, protocol=4)


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
            elif datatype == 'Ensembl_GTF':
                protein_id = seq.id.split('.')[0]
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

    def splitMutationByChromosome(self, chromosomes_genome_description=None, chromosomes_genome=None):
        '''split mutation file based on chromosomes
        '''
        tempfolder = self.tempfolder
        file_mutations = self.file_mutations
        if chromosomes_genome is None:
            chromosomes_genome = self.chromosomes_genome
        if chromosomes_genome_description is None:
            chromosomes_genome_description = self.chromosomes_genome_description

        df_mutations = pd.read_csv(file_mutations, sep='\t',low_memory=False)
        df_mutations['chr'] = df_mutations['chr'].astype(str)
        chromosomes_mutation = []
        
        for k,v in df_mutations.groupby('chr'):
            k_new = get_k_new(k, chromosomes_genome, chromosomes_genome_description)
            v = v.copy()
            v['chr'] = k_new
            # print(chromosomes_genome, k, k_new)
            tf = os.path.join(tempfolder, k_new + '.mutation.tsv')
            if k_new in chromosomes_genome:
                chromosomes_mutation.append(k_new)
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
            chromosome, feature, anno = es[0], es[2], es[-1]
            if feature not in features_keep:
                continue
            
            # parse anno
            tdc = parseGtfAttributes(anno)

            if protein_keyword in tdc:
                protein_ID = tdc[protein_keyword]
                if protein_ID in dc_protein2chr:
                    if dc_protein2chr[protein_ID] == chromosome:
                        # create output file if not exist
                        if chromosome not in dc_files:
                            tf = os.path.join(tempfolder, chromosome + '.gtf')
                            dc_files[chromosome] = open(tf,'w')
                        es[-1] = protein_ID
                        dc_files[chromosome].write('\t'.join(es) + '\n')
        
        for chromosome in dc_files:
            chromosomes_gtf.append(chromosome)
            dc_files[chromosome].close()
        print('finish splitting the gtf file')
        return chromosomes_gtf


    def getTranscript2ProteinRefSeq(self, dc_protein2chr):
        '''RefSeq, return dictionary of transcript to protein
        same proteins may exist in multiple chromosomes and each chromosome have different transcript_id, only keep the first chromosome in gtf file and the transcript_id in the first chromosome
        '''
        file_gtf = self.file_gtf
        protein_keyword = self.protein_keyword
        
        # check if already exist protein2chr
        file_protein2chr = os.path.join(self.tempfolder,'protein2chr.pickle')
        if os.path.exists(file_protein2chr):
            print(f'protein-chromosome is already determined in file {file_protein2chr}, use previous result')
            dc_protein2chr = pickle.load(open(file_protein2chr,'rb'))
            return dc_protein2chr
        
        # get transcript2protein
        dc_transcript2protein = {}

        for line in iterTxtWithComment(file_gtf):
            es = line.split('\t')
            chromosome, feature, anno = es[0], es[2], es[-1]
            if feature != 'CDS':
                continue
            tdc = parseGtfAttributes(anno)
            if protein_keyword in tdc:
                if 'transcript_id' in tdc:
                    if 'unknown_transcript' not in tdc['transcript_id']:
                        if dc_protein2chr[tdc[protein_keyword]] == chromosome:
                            dc_transcript2protein[tdc['transcript_id']] = tdc[protein_keyword]
                    else:
                        print('''RefSeq abnormal transcript_id "{}" for "{}"'''.format(tdc['transcript_id'], tdc[protein_keyword]))
        
        # write dc_protein2chr, as a checkpoint file
        with open(file_protein2chr,'wb') as f:
            pickle.dump(df, f, protocol=4)
        
        return dc_transcript2protein

    def splitGtfByChromosomesRefSeq(self,dc_protein2chr, dc_transcript2protein):
        '''split gtf file based on chromosome. only keep proteins in file_protein.
        For RefSeq, same proteins may exist in multiple chromosomes and each chromosome have different transcript_id, only keep the first chromosome in gtf file and the transcript_id in the first chromosome
        '''
        file_gtf = self.file_gtf
        tempfolder = self.tempfolder
        protein_keyword = self.protein_keyword

        chromosomes_gtf = []
        dc_files = {}
        features_keep = ['CDS', 'exon']
        
        for line in iterTxtWithComment(file_gtf):
            es = line.strip().split('\t')
            chromosome, feature, anno = es[0], es[2], es[-1]
            # filter based on feature
            if feature not in features_keep:
                continue

            tdc = parseGtfAttributes(anno)
            # for add protein_ID for exon rows
            if 'transcript_id' in tdc and protein_keyword not in tdc:
                transcript_id = tdc['transcript_id']
                if transcript_id in dc_transcript2protein:
                    tdc[protein_keyword] = dc_transcript2protein[transcript_id]
            
            if protein_keyword in tdc:
                protein_ID = tdc[protein_keyword]
                if protein_ID in dc_protein2chr:
                    if dc_protein2chr[protein_ID] == chromosome:
                        # create output file if not exist
                        if chromosome not in dc_files:
                            tf = os.path.join(tempfolder, chromosome + '.gtf')
                            dc_files[chromosome] = open(tf,'w')
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
        if self.file_mutations:
            self.chromosomes_mutation = self.splitMutationByChromosome()
        else:
            print('file_mutation is not provided, skip splitting mutation file')
        
        # split gtf
        if self.datatype != 'RefSeq':
            self.chromosomes_gtf = self.splitGtfByChromosomes(dc_protein2chr)
        else:
            dc_transcript2protein = self.getTranscript2ProteinRefSeq(dc_protein2chr)
            self.chromosomes_gtf = self.splitGtfByChromosomesRefSeq(dc_protein2chr, dc_transcript2protein)

        self.chromosomes = self.chromosomes_protein
        print('number of chromosomes with proteins:', len(self.chromosomes))
        return None
    
    def runSinglePerChrom(self, chromosome):
        '''run perChrom for single chromosome
        '''
        fchr_genome = os.path.join(self.tempfolder, chromosome + '.genome.fa')
        fchr_gtf = os.path.join(self.tempfolder, chromosome + '.gtf')
        fchr_mutations = os.path.join(self.tempfolder, chromosome + '.mutation.tsv')
        fchr_protein = os.path.join(self.tempfolder, chromosome + '.proteins.fa')

        # check if perChrom already run
        file_perChromFinished = os.path.join(self.tempfolder, chromosome + '.perChromFinished')
        if os.path.exists(file_perChromFinished):
            print('perChrom is already finished for chromosome', chromosome)
            return chromosome

        # check if all files exists
        if not os.path.exists(fchr_genome):
            print(chromosome, 'genome file not found. Proteins will be unchanged.')
            return None
        if not os.path.exists(fchr_gtf):
            print(chromosome, 'gtf file not found. Proteins will be unchanged.')
            return None
        if not os.path.exists(fchr_mutations):
            print(chromosome, 'mutation file not found. Proteins will be unchanged.')
            return None
        if not os.path.exists(fchr_protein):
            print(chromosome, 'protein file not found. Proteins will be unchanged.')
            return None
        
        # run PerChrom
        print('started running perChrom for chromosome:', chromosome)
        outprefix = os.path.join(self.tempfolder, chromosome)
        perchrom = PerChrom(file_genome = fchr_genome,
                    file_gtf = fchr_gtf,
                    file_mutations = fchr_mutations,
                    file_protein = fchr_protein,
                    threads = self.threads,
                    outprefix = outprefix,
                    datatype = self.datatype,
                    chromosome = chromosome)

        # run perChrom
        try:
            perchrom.run_perChrom()
            print('finished running perChrom for chromosome:', chromosome)
            open(file_perChromFinished,'w').write(str(chromosome))
            return chromosome
        except:
            print('cannot run perChrom for chromosome', chromosome, 'Proteins will be unchanged.')
            return None
        
        return None


    def runPerChom(self):
        '''run perChrom for all chromosomes
        '''
        chromosomes = self.chromosomes
        # run perChrom
        chromosomes_mutated = [self.runSinglePerChrom(e) for e in chromosomes]
        # successful chromosomes
        chromosomes_mutated = [e for e in chromosomes_mutated if e is not None]
        # failed chromosomes
        chromosomes_unchanged = [e for e in chromosomes if e not in chromosomes_mutated]
        # collect mutation annotations
        files_mutAnno = ['{}/{}.aa_mutations.csv'.format(self.tempfolder, chromosome) for chromosome in chromosomes_mutated]
        file_mutAnno = self.outprefix + '.pergeno.aa_mutations.csv'
        df_mutAnno = pd.concat([pd.read_csv(f, sep='\t') for f in files_mutAnno], ignore_index=True)
        print('total number of proteins with AA mutation:', df_mutAnno.shape[0])
        df_mutAnno.to_csv(file_mutAnno, sep='\t', index=None)

        # collect protein sequences
        files_proteins_changed = ['{}/{}.mutated_protein.fa'.format(self.tempfolder, chromosome) for chromosome in chromosomes_mutated]
        files_proteins_unchanged = ['{}/{}.proteins.fa'.format(self.tempfolder, chromosome) for chromosome in chromosomes_unchanged] + [os.path.join(self.tempfolder, 'not_in_genome.proteins.fa')]
        file_proteins_changed = self.outprefix + '.pergeno.protein_changed.fa'
        file_proteins_all = self.outprefix + '.pergeno.protein_all.fa'
        fout_proteins_changed = open(file_proteins_changed, 'w')
        fout_proteins_all = open(file_proteins_all, 'w')
        for f in files_proteins_changed:
            proteins_changed_ids = []
            for s in SeqIO.parse(f,'fasta'):
                fout_proteins_all.write('>{}\n{}\n'.format(s.description, str(s.seq)))
                if s.description.endswith('\tchanged'):
                    proteins_changed_ids.append(s.id)
                    fout_proteins_changed.write('>{}\n{}\n'.format(s.description, str(s.seq)))
            
            proteins_changed_ids = set(proteins_changed_ids)
            file_proteins_chromosome_input = f[:-19] + '.proteins.fa'
            for s in SeqIO.parse(file_proteins_chromosome_input,'fasta'):
                if s.id not in proteins_changed_ids:
                    fout_proteins_all.write('>{}\tunchanged\n{}\n'.format(s.description.split('\t',1)[1], str(s.seq)))

        
        for f in files_proteins_unchanged:
            for s in SeqIO.parse(f,'fasta'):
                fout_proteins_all.write('>{}\tunchanged\n{}\n'.format(s.description.split('\t',1)[1], str(s.seq)))

        fout_proteins_changed.close()
        fout_proteins_all.close()
        
        # clear temp folder
        if self.keep_all:
            print('keep all intermediate files')
        else:
            print('remove temp folder')
            shutil.rmtree(self.tempfolder)

        print('finished!')



description = '''PrecisionProDB_core, personal proteogenomic tools which outputs a new reference protein based on the variants data
'''

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-g','--genome', help = 'the reference genome sequence in fasta format. It can be a gzip file', required=True)
    parser.add_argument('-f', '--gtf', help='gtf file with CDS and exon annotations. It can be a gzip file', required=True)
    parser.add_argument('-m', '--mutations', help='a file stores the variants', required=True)
    parser.add_argument('-p','--protein', help = 'protein sequences in fasta format. It can be a gzip file. Only proteins in this file will be checked', required=True)
    parser.add_argument('-t', '--threads', help='number of threads/CPUs to run the program. default, use all CPUs available', type=int, default=os.cpu_count())
    parser.add_argument('-o', '--out', help='''output prefix. Three files will be saved, including the annotation for mutated transcripts, the mutated or all protein sequences. {out}.pergeno.aa_mutations.csv, {out}.pergeno.protein_all.fa, {out}.protein_changed.fa. default "perGeno" ''', default="perGeno")
    parser.add_argument('-a', '--datatype', help='''input datatype, could be GENCODE_GTF, GENCODE_GFF3, RefSeq, Ensembl_GTF or gtf. default "gtf". Ensembl_GFF3 is not supported. ''', default='gtf', type=str, choices=['GENCODE_GTF', 'GENCODE_GFF3','RefSeq','Ensembl_GTF','gtf'])
    parser.add_argument('-k','--protein_keyword', help='''field name in attribute column of gtf file to determine ids for proteins. default "auto", determine the protein_keyword based on datatype. "transcript_id" for GENCODE_GTF, "protein_id" for "RefSeq" and "Parent" for gtf and GENCODE_GFF3 ''', default='auto')
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
    keep_all = f.keep_all

    pergeno = PerGeno(
        file_genome = file_genome, 
        file_gtf=file_gtf, 
        file_mutations = file_mutations, 
        file_protein=file_protein, 
        threads=threads, 
        outprefix=outprefix, 
        datatype=datatype, 
        protein_keyword=protein_keyword, 
        keep_all = keep_all
        )
    # print(pergeno.__dict__)
    pergeno.splitInputByChromosomes()
    #print(pergeno.__dict__)
    pergeno.runPerChom()
