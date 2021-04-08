# Run PrecisonProDB for UniProt
`PATH_OF_PRECISONPRODB` is the location of PrecisionProDB.

Check this page for explanation of output files  
https://github.com/ATPs/PrecisionProDB/wiki/Outputs-of-PrecisionProDB

## Run with example files

### variants in tsv format
For UniProt models, PrecisonProDB will run for Ensembl models, and then do a match between UniProt proteins and Ensembl proteins to output the changed protein models.

```bash
cd $PATH_OF_PRECISONPRODB
cd examples

python ../src/PrecisionProDB.py -m gnomAD.variant.txt.gz -g ./Ensembl/Ensembl.genome.fa.gz -p ./Ensembl/Ensembl.protein.fa.gz -f ./Ensembl/Ensembl.gtf.gz -o ./UniProt/UniProt.tsv -a Ensembl_GTF --PEFF -U ./UniProt/UniProt.protein.fa.gz -D Uniprot

# compress files to save disk space
gzip ./UniProt/UniProt.tsv*
```

Here `--PEFF` means the PEFF option is enabled, and a result protein in [PEFF](http://www.psidev.info/peff) format with `VariantSimple` annotations will be generated.

Here, `-D` is set to be "Uniprot". When `-U`, `-m`, `-g`, `-p` and `-f` were all set, PrecisionProDB will skip the download step and use user-provided files. 

The output files are
* UniProt.tsv.pergeno.aa_mutations.csv.gz
* UniProt.tsv.pergeno.protein_all.fa.gz
* UniProt.tsv.pergeno.protein_changed.fa.gz
* UniProt.tsv.pergeno.protein_PEFF.fa.gz: the PEFF file.

Add additionally
* UniProt.tsv.uniprot_all.fa.gz
* UniProt.tsv.uniprot_changed.fa.gz
* UniProt.tsv.uniprot_changed.tsv.gz
* UniProt.tsv.uniprot_PEFF.fa.gz



### variants in vcf format
```bash
cd $PATH_OF_PRECISONPRODB
cd examples

python ../src/PrecisionProDB.py -m celline.vcf.gz -g ./Ensembl/Ensembl.genome.fa.gz -p ./Ensembl/Ensembl.protein.fa.gz -f ./Ensembl/Ensembl.gtf.gz -o ./UniProt/UniProt.vcf -a Ensembl_GTF --PEFF -U ./UniProt/UniProt.protein.fa.gz -D Uniprot

# compress files to save disk space
gzip ./UniProt/UniProt.vcf*
```

The output files are
* UniProt.vcf.pergeno.aa_mutations.csv.gz
* UniProt.vcf.pergeno.protein_all.fa.gz
* UniProt.vcf.pergeno.protein_changed.fa.gz
* UniProt.vcf.pergeno.protein_PEFF.fa.gz
* UniProt.vcf.vcf2mutation_1.tsv.gz
* UniProt.vcf.vcf2mutation_2.tsv.gz

Add additionally
* UniProt.vcf.uniprot_all.fa.gz
* UniProt.vcf.uniprot_changed.fa.gz
* UniProt.vcf.uniprot_changed.tsv.gz
* UniProt.vcf.uniprot_PEFF.fa.gz
