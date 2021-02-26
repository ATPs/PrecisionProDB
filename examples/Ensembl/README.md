# Run PrecisonProDB for Ensembl
`PATH_OF_PRECISONPRODB` is the location of PrecisionProDB.

Check this page for explanation of output files  
https://github.com/ATPs/PrecisionProDB/wiki/Outputs-of-PrecisionProDB

## Run with example files

### variants in tsv format
```bash
cd $PATH_OF_PRECISONPRODB
cd examples

python ../src/PrecisionProDB.py -m gnomAD.variant.txt.gz -g ./Ensembl/Ensembl.genome.fa.gz -p ./Ensembl/Ensembl.protein.fa.gz -f ./Ensembl/Ensembl.gtf.gz -o ./Ensembl/Ensembl.tsv -a Ensembl_GTF --PEFF

# compress files to save disk space
gzip ./Ensembl/Ensembl.tsv*
```

Here `--PEFF` means the PEFF option is enabled, and a result protein in [PEFF](http://www.psidev.info/peff) format with `VariantSimple` annotations will be generated.

The output files are
* Ensembl.tsv.pergeno.aa_mutations.csv.gz
* Ensembl.tsv.pergeno.protein_all.fa.gz
* Ensembl.tsv.pergeno.protein_changed.fa.gz
* Ensembl.tsv.pergeno.protein_PEFF.fa.gz: the PEFF file.


### variants in vcf format
```bash
cd $PATH_OF_PRECISONPRODB
cd examples

python ../src/PrecisionProDB.py -m celline.vcf.gz -g ./Ensembl/Ensembl.genome.fa.gz -p ./Ensembl/Ensembl.protein.fa.gz -f ./Ensembl/Ensembl.gtf.gz -o ./Ensembl/Ensembl.vcf -a Ensembl_GTF --PEFF

# compress files to save disk space
gzip ./Ensembl/Ensembl.vcf*
```

The output files are
* Ensembl.vcf.pergeno.aa_mutations.csv.gz
* Ensembl.vcf.pergeno.protein_all.fa.gz
* Ensembl.vcf.pergeno.protein_changed.fa.gz
* Ensembl.vcf.pergeno.protein_PEFF.fa.gz
* Ensembl.vcf.vcf2mutation_1.tsv.gz
* Ensembl.vcf.vcf2mutation_2.tsv.gz

