# Run PrecisonProDB for GENCODE
`PATH_OF_PRECISONPRODB` is the location of PrecisionProDB.

Check this page for explanation of output files  
https://github.com/ATPs/PrecisionProDB/wiki/Outputs-of-PrecisionProDB

## Run with example files

### variants in tsv format
```bash
cd $PATH_OF_PRECISONPRODB
cd examples

python ../src/PrecisionProDB.py -m gnomAD.variant.txt.gz -g ./GENCODE/GENCODE.genome.fa.gz -p ./GENCODE/GENCODE.protein.fa.gz -f ./GENCODE/GENCODE.gtf.gz -o ./GENCODE/GENCODE.tsv -a GENCODE_GTF --PEFF

# compress files to save disk space
gzip ./GENCODE/GENCODE.tsv*
```

Here `--PEFF` means the PEFF option is enabled, and a result protein in [PEFF](http://www.psidev.info/peff) format with `VariantSimple` annotations will be generated.

The output files are
* GENCODE.tsv.pergeno.aa_mutations.csv.gz
* GENCODE.tsv.pergeno.protein_all.fa.gz
* GENCODE.tsv.pergeno.protein_changed.fa.gz
* GENCODE.tsv.pergeno.protein_PEFF.fa.gz: the PEFF file.


### variants in vcf format
```bash
cd $PATH_OF_PRECISONPRODB
cd examples

python ../src/PrecisionProDB.py -m celline.vcf.gz -g ./GENCODE/GENCODE.genome.fa.gz -p ./GENCODE/GENCODE.protein.fa.gz -f ./GENCODE/GENCODE.gtf.gz -o ./GENCODE/GENCODE.vcf -a GENCODE_GTF --PEFF

# compress files to save disk space
gzip ./GENCODE/GENCODE.vcf*
```

The output files are
* GENCODE.vcf.pergeno.aa_mutations.csv.gz
* GENCODE.vcf.pergeno.protein_all.fa.gz
* GENCODE.vcf.pergeno.protein_changed.fa.gz
* GENCODE.vcf.pergeno.protein_PEFF.fa.gz
* GENCODE.vcf.vcf2mutation_1.tsv.gz
* GENCODE.vcf.vcf2mutation_2.tsv.gz

