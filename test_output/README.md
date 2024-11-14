# running the test scripts

file [testing.log](testing.log) contains the output of the test scripts and sample scripts of running PrecisonProDB in different modes.

file [test_running_summary.txt](test_running_summary.txt) is the summary of the test script outputs.

test scripts looks like

```bash
PATH_OF_PRECISIONPRODB=/data/p/xiaolong/PrecisionProDB
python $PATH_OF_PRECISIONPRODB/src/PrecisionProDB_test.py \
                -s $PATH_OF_PRECISIONPRODB \
                -o data/p/xiaolong/PrecisionProDB/test_output
```

then a combination of many testing runs were performed.

PrecisionProDB will be run with Ensembl, GENCODE, RefSeq, TransDecoder and UniProt gene annotations. 

It will use tsv, vcf or a variant string as input. 

It will run without building a sqlite file, building a sqlite file during the run or pre-building a sqlite file before running in the sqlite mode.

the outputs were shown in this folder.