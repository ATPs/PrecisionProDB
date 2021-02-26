# Example files to for testing PrecisionProDB

## Variant files
* a small vcf file: celline.vcf.gz
* a small tsv file: gnomAD.variant.txt.gz

For more larger files, users may find tsv files based on gnomAD 3.1 data from  
https://github.com/ATPs/PrecisionProDB_references/tree/main/gnomAD3.1

A bigger VCF file from Jurkat cellline can be found  
https://github.com/ATPs/PrecisionProDB_references/tree/main/cellline

## Run with only local files

Users can run PrecisionProDB with local genome files and gene models, which can is faster than the running mode of downloading the required files from the Internet.

Example files for each data types

* Ensembl
* GENCODE
* RefSeq
* TransDecoder
* UniProt

Check each of the subfolder to see how to run the testing codes.

## Run with the variant file only

PrecisonProDB can download the latest version of required files automatically from the Internet, so users can run get the results with the variant file only. However, as downloading large files may takesome time, the running time will be longer.

Check this page for explanation of output files  
https://github.com/ATPs/PrecisionProDB/wiki/Outputs-of-PrecisionProDB

### variants in tsv format
```
cd $PATH_OF_PRECISONPRODB
cd examples
python ../src/PrecisionProDB.py -m gnomAD.variant.txt.gz -D Ensembl -o ./PREFIX.tsv -a Ensembl_GTF --PEFF

```
Here `--PEFF` means the PEFF option is enabled, and a result protein in [PEFF](http://www.psidev.info/peff) format with `VariantSimple` annotations will be generated.

`-D` specified the gene models to use and to be downloaded. It can be `GENCODE`, `RefSeq`, `Ensembl` or `Uniprot`.

To run with your own data, just change `gnomAD.variant.txt.gz` to the location of your variant file.

Output files
* PREFIX.tsv.pergeno.aa_mutations.csv
* PREFIX.tsv.pergeno.protein_all.fa
* PREFIX.tsv.pergeno.protein_changed.fa
* PREFIX.tsv.pergeno.protein_PEFF.fa

### variants in vcf format

```
cd $PATH_OF_PRECISONPRODB
cd examples
python ../src/PrecisionProDB.py -m celline.vcf.gz -D Ensembl -o ./PREFIX.vcf -a Ensembl_GTF --PEFF

```
Output files:
* PREFIX.vcf.pergeno.aa_mutations.csv
* PREFIX.vcf.pergeno.protein_all.fa
* PREFIX.vcf.pergeno.protein_changed.fa
* PREFIX.vcf.pergeno.protein_PEFF.fa
* PREFIX.vcf.vcf2mutation_1.tsv
* PREFIX.vcf.vcf2mutation_2.tsv

Here, the `-m` is set to be `celline.vcf.gz`.

If the vcf provided by the user contails multiple samples, `-s/--sample` option can be used to specify a specific sample. If not set, the first sample in the vcf file will be used.

### Uniprot
For the two examples above, if `-D Uniprot` is set, several more output files will be generated.
* Uniprot.uniprot_all.fa
* Uniprot.uniprot_changed.fa
* Uniprot.uniprot_changed.tsv
* Uniprot.uniprot_PEFF.fa

These files should be used for downstream analysis.


