import os

file_genome = '/projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/GRCh38.p13.genome.fa.gz'
file_gtf = '/projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/gencode.v34.chr_patch_hapl_scaff.annotation.gtf.gz'
file_mutations = '/projectsp/f_jx76_1/xiaolong/2020humanRefPr/GENCODE/gnomAD3AF0.01ExonEthnics/adj.csv.gz'
file_proteins = '/projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/gencode.v34.pc_translations.fa.gz'
outprefix = '/projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/20200609gnomAD'

threads = os.cpu_count()
df_mutations = None
chromosomes = None
tempfolder = outprefix +'_temp'
datatype = 'gtf'
protein_keyword = 'transcript_id'

file_genome = '/projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/GRCh38.p13.genome.fa.gz'
file_gtf = '/projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/gencode.v34.chr_patch_hapl_scaff.annotation.gtf.gz'
file_mutations = '/projectsp/f_jx76_1/xiaolong/20200701perGeno/vcf/20200715_1KGP_selected.vcf.gz'
file_protein = '/projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/gencode.v34.pc_translations.fa.gz'
threads = 12
outprefix = '/projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.vcf'
datatype = 'gtf'
protein_keyword = 'auto'
filter_PASS = True
individual = None
chromosome_only = True

#python /cache/home/xc278/w/GitHub/perGeno/src/perGeno_core.py -g /projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/GRCh38.p13.genome.fa.gz -f /projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/gencode.v34.chr_patch_hapl_scaff.annotation.gtf.gz -m /projectsp/f_jx76_1/xiaolong/2020humanRefPr/GENCODE/gnomAD3AF0.01ExonEthnics/adj.csv.gz -p /projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/gencode.v34.pc_translations.fa.gz -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gtf 

#python /cache/home/xc278/w/GitHub/perGeno/src/perGeno_core.py -g /projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/GRCh38.p13.genome.fa.gz -f /projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/gencode.v34.chr_patch_hapl_scaff.annotation.gff3.gz -m /projectsp/f_jx76_1/xiaolong/2020humanRefPr/GENCODE/gnomAD3AF0.01ExonEthnics/adj.csv.gz -p /projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/gencode.v34.pc_translations.fa.gz -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gff3 

#python /cache/home/xc278/w/GitHub/perGeno/src/perGeno_core.py -g /projectsp/f_jx76_1/xiaolong/genome/human/RefSeq/GCF_000001405.39/GCF_000001405.39_GRCh38.p13_genomic.fna.gz -f /projectsp/f_jx76_1/xiaolong/genome/human/RefSeq/GCF_000001405.39/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz -m /projectsp/f_jx76_1/xiaolong/genome/gnomad/20200210AF/AFmostCommon/adj.csv.gz -p /projectsp/f_jx76_1/xiaolong/genome/human/RefSeq/GCF_000001405.39/GCF_000001405.39_GRCh38.p13_protein.faa.gz -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gtf

#python /cache/home/xc278/w/GitHub/perGeno/src/perGeno_core.py -g /projectsp/f_jx76_1/xiaolong/genome/human/RefSeq/GCF_000001405.39/GCF_000001405.39_GRCh38.p13_genomic.fna.gz -f /projectsp/f_jx76_1/xiaolong/genome/human/RefSeq/GCF_000001405.39/GCF_000001405.39_GRCh38.p13_genomic.gff.gz -m /projectsp/f_jx76_1/xiaolong/genome/gnomad/20200210AF/AFmostCommon/adj.csv.gz -p /projectsp/f_jx76_1/xiaolong/genome/human/RefSeq/GCF_000001405.39/GCF_000001405.39_GRCh38.p13_protein.faa.gz -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gff3

#python /cache/home/xc278/w/GitHub/perGeno/src/perGeno_core.py -g /projectsp/f_jx76_1/xiaolong/genome/human/CHESS/hg38_p8.fa -f /projectsp/f_jx76_1/xiaolong/20200701perGeno/CHESS/transdecoder/2020CHESS.StringTie.OneTranscriptPerGene.transcripts.translation.genome.gff3 -m /projectsp/f_jx76_1/xiaolong/genome/gnomad/20200210AF/AFmostCommon/adj.csv.gz -p /projectsp/f_jx76_1/xiaolong/20200701perGeno/CHESS/transdecoder/2020CHESS.StringTie.OneTranscriptPerGene.transcripts.fa.transdecoder.pep -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/CHESS/20200702CHESS.gff3 

#python /cache/home/xc278/w/GitHub/perGeno/src/perGeno_core.py -g /projectsp/f_jx76_1/xiaolong/genome/human/Ensembl/100/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -f /projectsp/f_jx76_1/xiaolong/genome/human/Ensembl/100/Homo_sapiens.GRCh38.100.gtf.gz -m /projectsp/f_jx76_1/xiaolong/2020humanRefPr/GENCODE/gnomAD3AF0.01ExonEthnics/adj.csv.gz -p /projectsp/f_jx76_1/xiaolong/genome/human/Ensembl/100/Homo_sapiens.GRCh38.pep.all.fa.gz -t 28 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/Ensembl//20200702Ensembl.gtf -a Ensembl_GTF



# python /cache/home/xc278/w/GitHub/perGeno/src/perGeno_chromosome.py -g /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gtf_temp/chr1.genome.fa -p /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gtf_temp/chr1.proteins.fa -f /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gtf_temp/chr1.gtf -m /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gtf_temp/chr1.mutation.tsv -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gtf_temp/chr1 -c chr1 -d GENCODE_GTF

# python /cache/home/xc278/w/GitHub/perGeno/src/perGeno_chromosome.py -g /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gff3_temp/chr1.genome.fa -p /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gff3_temp/chr1.proteins.fa -f /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gff3_temp/chr1.gtf -m /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gff3_temp/chr1.mutation.tsv -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gff3_temp/chr1 -c chr1 -d GENCODE_GFF3 

# python /cache/home/xc278/w/GitHub/perGeno/src/perGeno_chromosome.py -g /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gtf_temp/NC_000001.11.genome.fa -p /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gtf_temp/NC_000001.11.proteins.fa -f /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gtf_temp/NC_000001.11.gtf -m /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gtf_temp/NC_000001.11.mutation.tsv -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gtf_temp/NC_000001.11 -c NC_000001.11 -d RefSeq

# python /cache/home/xc278/w/GitHub/perGeno/src/perGeno_chromosome.py -g /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gff3_temp/NC_000001.11.genome.fa -p /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gff3_temp/NC_000001.11.proteins.fa -f /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gff3_temp/NC_000001.11.gtf -m /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gff3_temp/NC_000001.11.mutation.tsv -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gff3_temp/NC_000001.11 -c NC_000001.11 -d RefSeq

# python /cache/home/xc278/w/GitHub/perGeno/src/perGeno_chromosome.py -g /projectsp/f_jx76_1/xiaolong/20200701perGeno/CHESS/20200702CHESS.gff3_temp/chr1.genome.fa -p /projectsp/f_jx76_1/xiaolong/20200701perGeno/CHESS/20200702CHESS.gff3_temp/chr1.proteins.fa -f /projectsp/f_jx76_1/xiaolong/20200701perGeno/CHESS/20200702CHESS.gff3_temp/chr1.gtf -m /projectsp/f_jx76_1/xiaolong/20200701perGeno/CHESS/20200702CHESS.gff3_temp/chr1.mutation.tsv -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/CHESS/20200702CHESS.gff3_temp/chr1 -c chr1 -d gtf



# python /cache/home/xc278/w/GitHub/perGeno/src/perChrom.py -g /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gtf_temp/chr1.genome.fa -p /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gtf_temp/chr1.proteins.fa -f /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gtf_temp/chr1.gtf -m /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gtf_temp/chr1.mutation.tsv -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gtf_temp/chr1_2 -c chr1 -d GENCODE_GTF

# python /cache/home/xc278/w/GitHub/perGeno/src/perChrom.py -g /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gff3_temp/chr1.genome.fa -p /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gff3_temp/chr1.proteins.fa -f /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gff3_temp/chr1.gtf -m /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gff3_temp/chr1.mutation.tsv -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200702GENCODE.gff3_temp/chr1_2 -c chr1 -d GENCODE_GFF3 

# python /cache/home/xc278/w/GitHub/perGeno/src/perChrom.py -g /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gtf_temp/NC_000001.11.genome.fa -p /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gtf_temp/NC_000001.11.proteins.fa -f /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gtf_temp/NC_000001.11.gtf -m /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gtf_temp/NC_000001.11.mutation.tsv -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gtf_temp/NC_000001.11_2 -c NC_000001.11 -d RefSeq

# python /cache/home/xc278/w/GitHub/perGeno/src/perChrom.py -g /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gff3_temp/NC_000001.11.genome.fa -p /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gff3_temp/NC_000001.11.proteins.fa -f /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gff3_temp/NC_000001.11.gtf -m /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gff3_temp/NC_000001.11.mutation.tsv -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200702RefSeq.gff3_temp/NC_000001.11_2 -c NC_000001.11 -d RefSeq

# python /cache/home/xc278/w/GitHub/perGeno/src/perChrom.py -g /projectsp/f_jx76_1/xiaolong/20200701perGeno/CHESS/20200702CHESS.gff3_temp/chr1.genome.fa -p /projectsp/f_jx76_1/xiaolong/20200701perGeno/CHESS/20200702CHESS.gff3_temp/chr1.proteins.fa -f /projectsp/f_jx76_1/xiaolong/20200701perGeno/CHESS/20200702CHESS.gff3_temp/chr1.gtf -m /projectsp/f_jx76_1/xiaolong/20200701perGeno/CHESS/20200702CHESS.gff3_temp/chr1.mutation.tsv -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/CHESS/20200702CHESS.gff3_temp/chr1_2 -c chr1 -d gtf


#python /cache/home/xc278/w/GitHub/perGeno/src/perGeno_vcf.py -g /projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/GRCh38.p13.genome.fa.gz -f /projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/gencode.v34.chr_patch_hapl_scaff.annotation.gtf.gz -m /projectsp/f_jx76_1/xiaolong/20200701perGeno/vcf/20200715_1KGP_selected.vcf.gz -p /projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/gencode.v34.pc_translations.fa.gz -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200720GENCODE.gtf.vcf -s NA19704

#python /cache/home/xc278/w/GitHub/perGeno/src/perGeno_vcf.py -g /projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/GRCh38.p13.genome.fa.gz -f /projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/gencode.v34.chr_patch_hapl_scaff.annotation.gff3.gz -m /projectsp/f_jx76_1/xiaolong/20200701perGeno/vcf/20200715_1KGP_selected.vcf.gz -p /projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/gencode.v34.pc_translations.fa.gz -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/GENCODE/20200720GENCODE.gff3.vcf -s NA19704

#python /cache/home/xc278/w/GitHub/perGeno/src/perGeno_vcf.py -g /projectsp/f_jx76_1/xiaolong/genome/human/RefSeq/GCF_000001405.39/GCF_000001405.39_GRCh38.p13_genomic.fna.gz -f /projectsp/f_jx76_1/xiaolong/genome/human/RefSeq/GCF_000001405.39/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz -m /projectsp/f_jx76_1/xiaolong/20200701perGeno/vcf/20200715_1KGP_selected.vcf.gz -p /projectsp/f_jx76_1/xiaolong/genome/human/RefSeq/GCF_000001405.39/GCF_000001405.39_GRCh38.p13_protein.faa.gz -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200720RefSeq.gtf.vcf -s NA19704
#python /cache/home/xc278/w/GitHub/perGeno/src/perGeno_vcf.py -g /projectsp/f_jx76_1/xiaolong/genome/human/RefSeq/GCF_000001405.39/GCF_000001405.39_GRCh38.p13_genomic.fna.gz -f /projectsp/f_jx76_1/xiaolong/genome/human/RefSeq/GCF_000001405.39/GCF_000001405.39_GRCh38.p13_genomic.gff.gz -m /projectsp/f_jx76_1/xiaolong/20200701perGeno/vcf/20200715_1KGP_selected.vcf.gz -p /projectsp/f_jx76_1/xiaolong/genome/human/RefSeq/GCF_000001405.39/GCF_000001405.39_GRCh38.p13_protein.faa.gz -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/RefSeq/20200720RefSeq.gff3.vcf -s NA19704

#python /cache/home/xc278/w/GitHub/perGeno/src/perGeno_vcf.py -g /projectsp/f_jx76_1/xiaolong/genome/human/CHESS/hg38_p8.fa -f /projectsp/f_jx76_1/xiaolong/20200701perGeno/CHESS/transdecoder/2020CHESS.StringTie.OneTranscriptPerGene.transcripts.translation.genome.gff3 -m /projectsp/f_jx76_1/xiaolong/20200701perGeno/vcf/20200715_1KGP_selected.vcf.gz -p /projectsp/f_jx76_1/xiaolong/20200701perGeno/CHESS/transdecoder/2020CHESS.StringTie.OneTranscriptPerGene.transcripts.fa.transdecoder.pep -t 12 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/CHESS/20200720CHESS.gff3.vcf -s NA19704

#python /cache/home/xc278/w/GitHub/perGeno/src/perGeno_vcf.py -g /projectsp/f_jx76_1/xiaolong/genome/human/Ensembl/100/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -f /projectsp/f_jx76_1/xiaolong/genome/human/Ensembl/100/Homo_sapiens.GRCh38.100.gtf.gz -m /projectsp/f_jx76_1/xiaolong/20200701perGeno/vcf/20200715_1KGP_selected.vcf.gz -p /projectsp/f_jx76_1/xiaolong/genome/human/Ensembl/100/Homo_sapiens.GRCh38.pep.all.fa.gz -t 28 -o /projectsp/f_jx76_1/xiaolong/20200701perGeno/Ensembl//20200702Ensembl.vcf -a Ensembl_GTF 


chromosome = 'chr1'
file_genome = '/projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/20200609gnomAD_temp/chr1.genome.fa'
file_proteins = '/projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/20200609gnomAD_temp/chr1.proteins.fa'
file_gtf = '/projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/20200609gnomAD_temp/chr1.gtf'
file_mutations = '/projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/20200609gnomAD_temp/chr1.mutation.tsv'
outprefix = '/projectsp/f_jx76_1/xiaolong/genome/human/GENCODE/GRCh38.p13_release34/20200609gnomAD_temp/chr1'
cpu_counts = os.cpu_count()
datatype = 'GENCODE'


