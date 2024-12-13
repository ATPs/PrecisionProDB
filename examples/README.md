# Example files to for testing PrecisionProDB

## Variant files
* a small vcf file: celline.vcf.gz
* a small tsv file: gnomAD.variant.txt.gz

For more larger files, users may find tsv files based on gnomAD 3.1 data from  
https://github.com/ATPs/PrecisionProDB_references/tree/main/gnomAD3.1

A bigger VCF file from Jurkat cellline (including the one with T2T CHM13 as reference genome) can be found  
https://github.com/ATPs/PrecisionProDB_references/tree/main/cellline

[Jurkat.CHM13.RefSeq.vcf.gz](https://github.com/ATPs/PrecisionProDB_references/blob/main/cellline/Jurkat.CHM13.RefSeq.vcf.gz)

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

PrecisonProDB can download the latest version of required files automatically from the Internet, so users can run get the results with the variant file only. However, as downloading large files may take some time, the running time will be longer.

Check this page for explanation of output files  
https://github.com/ATPs/PrecisionProDB/wiki/Outputs-of-PrecisionProDB

### enable sqlite
```bash
cd $PATH_OF_PRECISONPRODB
cd examples
python ../src/PrecisionProDB.py -m gnomAD.variant.txt.gz -D Ensembl -o ./PREFIX.tsv -a Ensembl_GTF --PEFF \
                    -S path_to_file_sqlite

```
If `path_to_file_sqlite` does not exist, a sqlite file will be created.

If `path_to_file_sqlite` is provided, the program will use informations stored in the sqlite file to speed up.

file `path_to_file_sqlite` can be build in advance with the following command.

```bash
cd $PATH_OF_PRECISONPRODB
python ./src/buildSqlite.py  -g GRCh38.p14.genome.fa.gz \
                    -p gencode.v47.pc_translations.fa.gz \
                    -f gencode.v47.chr_patch_hapl_scaff.annotation.gtf.gz \
                    -o GENCODE.tsv.buildSqlite \
                    -a GENCODE_GTF \
                    -S GENCODE.sqlite \
                    -t 8
```

### variants in tsv format
```bash
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

### variant in string format

```Python
cd $PATH_OF_PRECISONPRODB

python ./src/PrecisionProDB.py \
        -m chr1-942451-T-C,1-6253878-C-T,1-2194700-C-G,1-1719406-G-A \
        -o ./test_output/Ensembl/str/sqlite_two_step/Ensembl.str.sqlite_two_step \
        -a Ensembl_GTF --PEFF -t 4 \
        -S ./test_output/Ensembl/Ensembl.sqlite
```
two output will be generated

* Ensembl.str.sqlite_two_step.pergeno.aa_mutations.csv
* Ensembl.str.sqlite_two_step.pergeno.mutated_protein.fa

Note:
if variant is a string, the `-S` option is required.

## CHM13 T2T from RefSeq

Since the data is in the same format like RefSeq, we tested the command below and it works perfectly.

```Python
cd $PATH_OF_PRECISONPRODB

python ./src/downloadHuman.py -d CHM13 -o precisionprodb

python ./src/PrecisionProDB.py \
        -g precisionprodb/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz \
        -p precisionprodb/GCF_009914755.1_T2T-CHM13v2.0_protein.faa.gz \
        -t 12 \
        -a RefSeq \
        -f precisionprodb/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz \
        -m precisionprodb/Jurkat.CHM13.RefSeq.vcf.gz \
        -o precisionprodb/Jurkat.CHM13.RefSeq.PrecisonProDB \
        -S precisionprodb/GCF_009914755.1_T2T-CHM13v2.0.sqlite
```

# vcf2mutation
to get tsv file from the vcf file, you can use the vcf2mutation.py script
```
cd $PATH_OF_PRECISONPRODB
python ./src/vcf2mutation.py -h

usage: vcf2mutation [-h] -i FILE_VCF [-o OUTPREFIX] [-s SAMPLE]
                    [-F] [-A]

convert extract mutation information from vcf file

options:
  -h, --help            show this help message and exit
  -i FILE_VCF, --file_vcf FILE_VCF
                        input vcf file. It can be a gzip file
  -o OUTPREFIX, --outprefix OUTPREFIX
                        output prefix to store the two output
                        dataframes, default: None, do not write
                        the result to files. file will be
                        outprefix_1/2.tsv
  -s SAMPLE, --sample SAMPLE
                        sample name in the vcf to extract the
                        variant information. default: None,
                        extract the first sample
  -F, --no_filter       default only keep variant with value
                        "PASS" FILTER column of vcf file. if
                        set, do not filter
  -A, --all_chromosomes
                        default keep variant in chromosomes and
                        ignore those in short fragments of the
                        genome. if set, use all chromosomes
                        including fragments

```
