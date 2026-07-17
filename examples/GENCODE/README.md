# Run PrecisonProDB for GENCODE
`PATH_OF_PRECISONPRODB` is the location of PrecisionProDB.

Check this page for explanation of output files  
https://github.com/ATPs/PrecisionProDB/wiki/Outputs-of-PrecisionProDB

## Run with example files

### variants in tsv format
```bash
cd $PATH_OF_PRECISONPRODB
cd examples

python ../src/precisionprodb/PrecisionProDB.py -m gnomAD.variant.txt.gz -g ./GENCODE/GENCODE.genome.fa.gz -p ./GENCODE/GENCODE.protein.fa.gz -f ./GENCODE/GENCODE.gtf.gz -o ./GENCODE/GENCODE.tsv -a GENCODE_GTF --PEFF

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

python ../src/precisionprodb/PrecisionProDB.py -m celline.vcf.gz -g ./GENCODE/GENCODE.genome.fa.gz -p ./GENCODE/GENCODE.protein.fa.gz -f ./GENCODE/GENCODE.gtf.gz -o ./GENCODE/GENCODE.vcf -a GENCODE_GTF --PEFF

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

## Novel peptide output

`--peptide` reports altered-protein digestion products that are absent from the canonical GENCODE proteome under the same digestion profile. It requires a real annotation SQLite file; do not use `--sqlite NONE`.

Run the following commands from the PrecisionProDB repository root. The first command builds the annotation SQLite. The second command enables peptide mode and automatically creates a profile-specific canonical peptide SQLite index beside it when one is not already present.

```bash
python src/precisionprodb/buildSqlite.py \
  -S examples/GENCODE/GENCODE.peptide.sqlite \
  -g examples/GENCODE/GENCODE.genome.fa.gz \
  -p examples/GENCODE/GENCODE.protein.fa.gz \
  -f examples/GENCODE/GENCODE.gtf.gz \
  -a GENCODE_GTF \
  -o examples/GENCODE/GENCODE.peptide.build \
  -t 2

python src/precisionprodb/PrecisionProDB.py \
  -m examples/celline.vcf.gz \
  -S examples/GENCODE/GENCODE.peptide.sqlite \
  -o examples/GENCODE/GENCODE.peptide \
  -a GENCODE_GTF \
  --peptide \
  -t 2

sed -n '1,4p' examples/GENCODE/GENCODE.peptide.pergeno.peptide_novel.tsv
sed -n '1,8p' examples/GENCODE/GENCODE.peptide.pergeno.peptide_novel.fa
```

In addition to the ordinary VCF outputs, this creates the following uncompressed files:

* `GENCODE.peptide.pergeno.peptide_novel.tsv`: peptide-to-altered-protein mappings with mutation and digestion-profile annotations.
* `GENCODE.peptide.pergeno.peptide_novel.fa`: nonredundant peptide sequences with `PPDBpep_######|n_mappings=N` headers.

The default profile uses Trypsin, two missed cleavages, lengths 7–35, full specificity, both initiator-methionine forms, and exact I/L comparison. See the [Novel peptide output](https://github.com/ATPs/PrecisionProDB/wiki/Novel-peptide-output) Wiki page for custom profiles, reusable indexes, and complete output-column descriptions.
