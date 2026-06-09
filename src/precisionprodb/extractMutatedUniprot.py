from Bio import SeqIO
import pandas as pd

# files_uniprot = '/projectsp/f_jx76_1/xiaolong/2020humanRefPr/Uniprot/20201212Test/UP000005640_9606.fasta.gz,/projectsp/f_jx76_1/xiaolong/2020humanRefPr/Uniprot/20201212Test/UP000005640_9606_additional.fasta.gz'
# files_ref = '/projectsp/f_jx76_1/xiaolong/genome/human/Ensembl/101/Homo_sapiens.GRCh38.pep.all.fa.gz,/projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release35/gencode.v35.pc_translations.fa.gz'
# files_alt = '/projectsp/f_jx76_1/xiaolong/2020humanRefPr/Ensembl/20201119EthnicProteins/Ensembl_adj.pergeno.protein_all.fa,/projectsp/f_jx76_1/xiaolong/2020humanRefPr/GENCODE/20201119EthnicProteins/GENCODE_adj.pergeno.protein_all.fa'
# outprefix = '/projectsp/f_jx76_1/xiaolong/2020humanRefPr/Uniprot/20201212Test/20201212TestUniprot.perGeno'
# length_min = 20

def openFile(filename):
    '''open text or gzipped text file
    '''
    if filename.endswith('.gz'):
        import gzip
        return gzip.open(filename,'rt')
    return open(filename)

def extractMutatedUniprot(files_uniprot, files_ref, files_alt, outprefix, length_min = 20):
    '''match uniprot proteins in files_uniprot with files_ref. To start a match, the min length of protein is length_min.
    output mutated uniprot proteins in files_alt.
    write three files, outprefix + '.uniprot_changed.tsv'/'.uniprot_changed.fa'/'.uniprot_all.fa'
    '''

    df = pd.DataFrame()
    # store uniprot sequences to a dataframe
    ls_files = files_uniprot.split(',')
    ls_seqs = []
    for f in ls_files:
        for seq in SeqIO.parse(openFile(f),'fasta'):
            ls_seqs.append(seq)
    df['uniprot_id'] = [e.id for e in ls_seqs]
    df['uniprot_description'] = [e.description for e in ls_seqs]
    df['uniprot_seq_ori'] = [str(e.seq).strip('*X').upper() for e in ls_seqs]

    # get reference protein_ids
    dc_seq2ref = dict.fromkeys([e for e in list(df['uniprot_seq_ori']) if len(e) >= length_min])
    ls_files = files_ref.split(',')
    for f in ls_files:
        for s in SeqIO.parse(openFile(f),'fasta'):
            protein_seq = str(s.seq).strip('*X').upper()
            if protein_seq in dc_seq2ref:
                dc_seq2ref[protein_seq] = s.id

    df['ref_id'] = df['uniprot_seq_ori'].map(dc_seq2ref)
    # filter to include only uniprot_id with ref_id
    df1 = df[df['ref_id'].notnull()]

    # for ref_id in df1, only store those with mutations based on files_alt
    dc_ref2alt = {k:set() for k in df1['ref_id']}
    ls_files = files_alt.split(',')
    for f in ls_files:
        for s in SeqIO.parse(openFile(f),'fasta'):
            if s.description.endswith('\tchanged'):
                protein_id = s.id.split('__')[0]
                if protein_id in dc_ref2alt:
                    protein_seq = str(s.seq).strip('*X').upper()
                    dc_ref2alt[protein_id].add(protein_seq)
    dc_ref2alt = {k:v for k,v in dc_ref2alt.items() if len(v) > 0}

    df2 = df1[df1['ref_id'].isin(dc_ref2alt)]

    # save relationship of uniprot2ref
    df3 = df2[['uniprot_id', 'ref_id']]
    df3.to_csv(outprefix + '.uniprot_changed.tsv', sep='\t', index=None)

    # save changed proteins and all proteins
    file_changed = outprefix + '.uniprot_changed.fa'
    file_all = outprefix + '.uniprot_all.fa'
    fout_changed = open(file_changed,'w')
    fout_all = open(file_all, 'w')

    for _, r in df2.iterrows():
        uniprot_id, uniprot_description, ref_id = r['uniprot_id'], r['uniprot_description'], r['ref_id']
        alt_seqs = list(dc_ref2alt[ref_id])
        uniprot_description = uniprot_description[len(uniprot_id):].strip()
        if len(alt_seqs) == 1:
            fout_changed.write('>{} {}\tchanged\n{}\n'.format(uniprot_id,uniprot_description,alt_seqs[0]))
            fout_all.write('>{} {}\tchanged\n{}\n'.format(uniprot_id,uniprot_description,alt_seqs[0]))
        else:
            for n, alt_seq in enumerate(alt_seqs):
                fout_changed.write('>{}__{} {}\tchanged\n{}\n'.format(uniprot_id,n+1, uniprot_description,alt_seq))
                fout_all.write('>{}__{} {}\tchanged\n{}\n'.format(uniprot_id,n+1, uniprot_description,alt_seq))

    changed_ids = set(df2['uniprot_id'])
    for _, r in df.iterrows():
        uniprot_id = r['uniprot_id']
        if uniprot_id not in changed_ids:
            uniprot_description = r['uniprot_description']
            uniprot_seq = r['uniprot_seq_ori']
            fout_all.write('>{}\tunchanged\n{}\n'.format(uniprot_description, uniprot_seq))

    fout_changed.close()
    fout_all.close()




description = '''
Match UniProt proteins in files_uniprot with files_ref.
Output mutated proteins in files_alt.
write three files, outprefix + '.uniprot_changed.tsv'/'.uniprot_changed.fa'/'.uniprot_all.fa'
'''

def main():
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-u','--files_uniprot', help = 'uniprot proteins. If more than one files, join by ","', required=True)
    parser.add_argument('-r', '--files_ref', help='reference proteins to match with uniprot proteins. If more than one files, join by ","')
    parser.add_argument('-a', '--files_alt', help='altered reference proteins. If more than one files, join by ",". The order should be the same as files_ref')
    parser.add_argument('-o', '--outprefix', help='prefix for output files. default:"perGeno"', default='perGeno')
    parser.add_argument('-m', '--length_min', help='minumum length required when matching UniProt sequences with sequences in files_ref. default: "20"', default=20, type=int)
    f = parser.parse_args()
    
    files_uniprot = f.files_uniprot
    files_ref = f.files_ref
    files_alt = f.files_alt
    outprefix = f.outprefix
    length_min = f.length_min
    extractMutatedUniprot(files_uniprot, files_ref, files_alt, outprefix, length_min = length_min)
if __name__ == '__main__':
    main()