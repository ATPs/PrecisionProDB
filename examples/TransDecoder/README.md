# Run PrecisonProDB for TransDecoder

## Translation of gene models with TransDecoder

Check the following link for details.  
https://github.com/ATPs/PrecisionProDB/wiki/An-Example-of-Running-TransDecoder

Running TransDecoder with example files

We provided example files for running TransDecoder. The files were stored in `PATH_OF_PRECISONPRODB/examples/TransDecoder`.

`PATH_OF_PRECISONPRODB` is the location of PrecisionProDB.  
`PATH_OF_TRANSDECODER` is the location of TransDecoder.

```bash
# change working directory to example folder
cd PATH_OF_PRECISONPRODB/examples/TransDecoder

# decompress the example files
gzip -d *.gz

# get transcript sequences from the gff file
$PATH_OF_TRANSDECODER/util/gtf_genome_to_cdna_fasta.pl TransDecoder.gtf TransDecoder.genome.fa >TransDecoder.transcripts.fa

# convert gff file to gff3 format
$PATH_OF_TRANSDECODER/util/gtf_to_alignment_gff3.pl TransDecoder.gtf >TransDecoder.gff3

# translate and predict proteins. -m is the minimum length of proteins.
$PATH_OF_TRANSDECODER/TransDecoder.LongOrfs -t TransDecoder.transcripts.fa -m 60
$PATH_OF_TRANSDECODER/TransDecoder.Predict -t TransDecoder.transcripts.fa

# map translated proteins to the genome
$PATH_OF_TRANSDECODER/util/cdna_alignment_orf_to_genome_orf.pl \
                            TransDecoder.transcripts.fa.transdecoder.gff3 \
                            TransDecoder.gff3 \
                            TransDecoder.transcripts.fa > TransDecoder.transcripts.fa.transdecoder.genome.gff3

# to save disk space, compress these files
gzip *
```

The output files are
* TransDecoder.gff3
* TransDecoder.transcripts.fa
* TransDecoder.transcripts.fa.transdecoder.bed
* TransDecoder.transcripts.fa.transdecoder.cds
* TransDecoder.transcripts.fa.transdecoder.genome.gff3: gff3 annotation based on the genome sequence with proteins.
* TransDecoder.transcripts.fa.transdecoder.gff3: gff3 annotation of proteins based on transcript sequences.
* TransDecoder.transcripts.fa.transdecoder.pep: final result of translated proteins.

Of these files, `TransDecoder.transcripts.fa.transdecoder.pep` and `TransDecoder.transcripts.fa.transdecoder.genome.gff3` are inputs for PrecisonProDB.

You might see some warning messages, which are just normal warning messages.

`PATH_OF_PRECISONPRODB` is the location of PrecisionProDB.

Check this page for explanation of output files  
https://github.com/ATPs/PrecisionProDB/wiki/Outputs-of-PrecisionProDB

## Run with example files

### variants in tsv format
```bash
cd $PATH_OF_PRECISONPRODB
cd examples

python ../src/PrecisionProDB.py -m gnomAD.variant.txt.gz -g ./TransDecoder/TransDecoder.genome.fa.gz -p ./TransDecoder/TransDecoder.transcripts.fa.transdecoder.pep.gz -f ./TransDecoder/TransDecoder.transcripts.fa.transdecoder.genome.gff3.gz -o ./TransDecoder/TransDecoder.tsv -a gtf --PEFF

# compress files to save disk space
gzip ./TransDecoder/TransDecoder.tsv*
```

Here `--PEFF` means the PEFF option is enabled, and a result protein in [PEFF](http://www.psidev.info/peff) format with `VariantSimple` annotations will be generated.

The output files are
* TransDecoder.tsv.pergeno.aa_mutations.csv.gz
* TransDecoder.tsv.pergeno.protein_all.fa.gz
* TransDecoder.tsv.pergeno.protein_changed.fa.gz
* TransDecoder.tsv.pergeno.protein_PEFF.fa.gz: the PEFF file.


### variants in vcf format
```bash
cd $PATH_OF_PRECISONPRODB
cd examples

python ../src/PrecisionProDB.py -m celline.vcf.gz -g ./TransDecoder/TransDecoder.genome.fa.gz -p ./TransDecoder/TransDecoder.transcripts.fa.transdecoder.pep.gz -f ./TransDecoder/TransDecoder.transcripts.fa.transdecoder.genome.gff3.gz -o ./TransDecoder/TransDecoder.vcf -a gtf --PEFF

# compress files to save disk space
gzip ./TransDecoder/TransDecoder.vcf*
```

The output files are
* TransDecoder.vcf.pergeno.aa_mutations.csv.gz
* TransDecoder.vcf.pergeno.protein_all.fa.gz
* TransDecoder.vcf.pergeno.protein_changed.fa.gz
* TransDecoder.vcf.pergeno.protein_PEFF.fa.gz
* TransDecoder.vcf.vcf2mutation_1.tsv.gz
* TransDecoder.vcf.vcf2mutation_2.tsv.gz

