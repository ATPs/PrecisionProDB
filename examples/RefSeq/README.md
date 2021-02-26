# Run PrecisonProDB for RefSeq
`PATH_OF_PRECISONPRODB` is the location of PrecisionProDB.

Check this page for explanation of output files  
https://github.com/ATPs/PrecisionProDB/wiki/Outputs-of-PrecisionProDB

## Run with example files

### variants in tsv format
```bash
cd $PATH_OF_PRECISONPRODB
cd examples

python ../src/PrecisionProDB.py -m gnomAD.variant.txt.gz -g ./RefSeq/RefSeq.genome.fa.gz -p ./RefSeq/RefSeq.protein.fa.gz -f ./RefSeq/RefSeq.gtf.gz -o ./RefSeq/RefSeq.tsv -a RefSeq --PEFF

# compress files to save disk space
gzip ./RefSeq/RefSeq.tsv*
```

Here `--PEFF` means the PEFF option is enabled, and a result protein in [PEFF](http://www.psidev.info/peff) format with `VariantSimple` annotations will be generated.

The output files are
* RefSeq.tsv.pergeno.aa_mutations.csv.gz
* RefSeq.tsv.pergeno.protein_all.fa.gz
* RefSeq.tsv.pergeno.protein_changed.fa.gz
* RefSeq.tsv.pergeno.protein_PEFF.fa.gz: the PEFF file.


### variants in vcf format
```bash
cd $PATH_OF_PRECISONPRODB
cd examples

python ../src/PrecisionProDB.py -m celline.vcf.gz -g ./RefSeq/RefSeq.genome.fa.gz -p ./RefSeq/RefSeq.protein.fa.gz -f ./RefSeq/RefSeq.gtf.gz -o ./RefSeq/RefSeq.vcf -a RefSeq --PEFF

# compress files to save disk space
gzip ./RefSeq/RefSeq.vcf*
```

The output files are
* RefSeq.vcf.pergeno.aa_mutations.csv.gz
* RefSeq.vcf.pergeno.protein_all.fa.gz
* RefSeq.vcf.pergeno.protein_changed.fa.gz
* RefSeq.vcf.pergeno.protein_PEFF.fa.gz
* RefSeq.vcf.vcf2mutation_1.tsv.gz
* RefSeq.vcf.vcf2mutation_2.tsv.gz

