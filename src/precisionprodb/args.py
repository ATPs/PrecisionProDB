import os


DEFAULT_THREADS = min(20, os.cpu_count() or 1)
DATATYPE_CHOICES = ['GENCODE_GTF', 'GENCODE_GFF3', 'RefSeq', 'Ensembl_GTF', 'gtf']
DOWNLOAD_CHOICES = ['GENCODE', 'RefSeq', 'Ensembl', 'Uniprot', 'CHM13', '']


ARGUMENT_SPECS = {
    'genome': {
        'flags': ('-g', '--genome'),
        'kwargs': {
            'help': 'the reference genome sequence in fasta format. It can be a gzip file',
            'default': '',
        },
    },
    'gtf': {
        'flags': ('-f', '--gtf'),
        'kwargs': {
            'help': 'gtf file with CDS and exon annotations. It can be a gzip file',
            'default': '',
        },
    },
    'protein': {
        'flags': ('-p', '--protein'),
        'kwargs': {
            'help': 'protein sequences in fasta format. It can be a gzip file. Only proteins in this file will be checked',
            'default': '',
        },
    },
    'mutations_detailed': {
        'flags': ('-m', '--mutations'),
        'kwargs': {
            'help': '''
                        a file stores the variants.
                        If the file ends with ".vcf" or ".vcf.gz", treat as vcf input. Otherwise, treat as TSV input.
                        A string like "chr1-788418-CAG-C" or "chr1-942451-T-C,1-6253878-C-T,1-2194700-C-G" can used as variant input, too. In this mode, --sample will not be used.
                        If multiple vcf files are provided, use "," to join the file names. For example, "--mutations file1.vcf,file2.vcf". A pattern match is also supported for input vcf, but quote is required to get it work. For example '--mutations "file*.vcf" '
                        A TSV manifest with first header column "filepath" is supported for one-VCF-per-sample population mode. Optional manifest columns: sample, name_use.

                        ''',
            'default': '',
            'required': False,
        },
    },
    'mutations_basic': {
        'flags': ('-m', '--mutations'),
        'kwargs': {
            'help': 'a file stores the variants',
            'required': True,
        },
    },
    'mutations_perchrom': {
        'flags': ('-m', '--mutations'),
        'kwargs': {
            'help': 'a file stores the variants, or a string looks like "1-55051215-G-GA" or "chr1-55051215-G-GA" ',
            'required': True,
        },
    },
    'threads': {
        'flags': ('-t', '--threads'),
        'kwargs': {
            'help': 'number of threads/CPUs to run the program. default, use 20 or all CPUs available, whichever is smaller',
            'type': int,
            'default': DEFAULT_THREADS,
        },
    },
    'threads_vcf2mutation': {
        'flags': ('-t', '--threads'),
        'kwargs': {
            'help': 'number of threads/CPUs to run the program. default, 1. For manifest input, VCF files are parsed in parallel across manifest rows.',
            'type': int,
            'default': 1,
        },
    },
    'out': {
        'flags': ('-o', '--out'),
        'kwargs': {
            'help': '''output prefix, folder path could be included. Three or five files will be saved depending on the variant file format. Outputs include the annotation for mutated transcripts, the mutated or all protein sequences, two variant files from vcf. {out}.pergeno.aa_mutations.csv, {out}.pergeno.protein_all.fa, {out}.protein_changed.fa, {out}.vcf2mutation_1/2.tsv. default "perGeno" ''',
            'default': 'perGeno',
        },
    },
    'datatype': {
        'flags': ('-a', '--datatype'),
        'kwargs': {
            'help': '''input datatype, could be GENCODE_GTF, GENCODE_GFF3, RefSeq, Ensembl_GTF or gtf. default "gtf". Ensembl_GFF3 is not supported. ''',
            'default': 'gtf',
            'type': str,
            'choices': DATATYPE_CHOICES,
        },
    },
    'protein_keyword': {
        'flags': ('-k', '--protein_keyword'),
        'kwargs': {
            'help': '''field name in attribute column of gtf file to determine ids for proteins. default "auto", determine the protein_keyword based on datatype. "transcript_id" for GENCODE_GTF, "protein_id" for "RefSeq" and "Parent" for gtf and GENCODE_GFF3 ''',
            'default': 'auto',
        },
    },
    'no_filter': {
        'flags': ('-F', '--no_filter'),
        'kwargs': {
            'help': 'default only keep variant with value "PASS" FILTER column of vcf file. if set, do not filter',
            'action': 'store_true',
        },
    },
    'sample_mixed': {
        'flags': ('-s', '--sample'),
        'kwargs': {
            'help': '''
                        sample name in the vcf/tsv to extract the variant information. default: None, extract the first sample in vcf file, or use all variants in the tsv file.
                        For multiple samples, use "," to join the sample names. For example, "--sample sample1,sample2,sample3".
                        To use all samples, use "--sample ALL_SAMPLES".
                        To use all variants regardless where the variants from, use "--sample ALL_VARIANTS".
                        ''',
            'default': None,
        },
    },
    'sample_vcf': {
        'flags': ('-s', '--sample'),
        'kwargs': {
            'help': 'sample name in the vcf to extract the variant information. default: None, extract the first sample. For multiple samples, use "," to join the sample names. For example, "--sample sample1,sample2,sample3". To use all samples, use "--sample ALL_SAMPLES". To use all variants regardless where the variants from, use "--sample ALL_VARIANTS".',
            'default': None,
        },
    },
    'sample_perchrom': {
        'flags': ('-s', '--sample'),
        'kwargs': {
            'help': '''
                        sample name in the vcf to extract the variant information. default: None, use all variants and do not consider samples.
                        For multiple samples, use "," to join the sample names. For example, "--sample sample1,sample2,sample3".
                        To use all samples, use "--sample ALL_SAMPLES".
                        To use all variants regardless where the variants from, use "--sample ALL_VARIANTS".
                        ''',
            'default': None,
        },
    },
    'all_chromosomes': {
        'flags': ('-A', '--all_chromosomes'),
        'kwargs': {
            'help': 'default keep variant in chromosomes and ignore those in short fragments of the genome. if set, use all chromosomes including fragments when parsing the vcf file',
            'action': 'store_true',
        },
    },
    'download': {
        'flags': ('-D', '--download'),
        'kwargs': {
            'help': '''download could be 'GENCODE','RefSeq','Ensembl','Uniprot', 'CHM13'. If set, PrecisonProDB will try to download genome, gtf and protein files from the Internet. Download will be skipped if "--genome, --gtf, --protein, (--uniprot)" were all set. Settings from "--genome, --gtf, --protein, (--uniprot), --datatype" will not be used if the files were downloaded by PrecisonProDB. default "". Note, if --sqlite is set, will not download any files ''',
            'default': '',
            'type': str,
            'choices': DOWNLOAD_CHOICES,
        },
    },
    'uniprot': {
        'flags': ('-U', '--uniprot'),
        'kwargs': {
            'help': '''uniprot protein sequences. If more than one file, use "," to join the files. default "". For example, "UP000005640_9606.fasta.gz", or "UP000005640_9606.fasta.gz,UP000005640_9606_additional.fasta" ''',
            'default': '',
            'type': str,
        },
    },
    'uniprot_min_len': {
        'flags': ('--uniprot_min_len',),
        'kwargs': {
            'help': '''minimum length required when matching uniprot sequences to proteins annotated in the genome. default 20 ''',
            'default': 20,
            'type': int,
        },
    },
    'peff': {
        'flags': ('--PEFF',),
        'kwargs': {
            'help': 'If set, PEFF format file(s) will be generated. Default: do not generate PEFF file(s).',
            'action': 'store_true',
        },
    },
    'keep_all': {
        'flags': ('--keep_all',),
        'kwargs': {
            'help': 'If set, do not delete files generated during the run',
            'action': 'store_true',
        },
    },
    'sqlite': {
        'flags': ('-S', '--sqlite'),
        'kwargs': {
            'help': '''A path of sqlite file for re-use of annotation info. default outprefix + '.sqlite'. The program will create a sqlite file if the file does not exist. If the file already exists, the program will use data stored in the file. It will cause error if the content in the sqlite file is not as expected. To disable sqlite, set to "NONE". ''',
            'default': '',
            'type': str,
        },
    },
    'peptide': {
        'flags': ('--peptide',),
        'kwargs': {
            'help': 'Generate globally novel altered-protein peptides using a canonical peptide SQLite index.',
            'action': 'store_true',
        },
    },
    'peptide_sqlite': {
        'flags': ('--peptide-sqlite',),
        'kwargs': {
            'help': 'Canonical peptide SQLite index. Default: a profile-specific path beside --sqlite.',
            'default': '',
        },
    },
    'peptide_enzyme': {
        'flags': ('--peptide-enzyme',),
        'kwargs': {'default': 'Trypsin', 'help': 'Peptide digestion enzyme. Default: Trypsin.'},
    },
    'peptide_missed_cleavages': {
        'flags': ('--peptide-missed-cleavages',),
        'kwargs': {'type': int, 'default': 2, 'help': 'Maximum missed cleavages. Default: 2.'},
    },
    'peptide_min_length': {
        'flags': ('--peptide-min-length',),
        'kwargs': {'type': int, 'default': 7, 'help': 'Minimum peptide length. Default: 7.'},
    },
    'peptide_max_length': {
        'flags': ('--peptide-max-length',),
        'kwargs': {'type': int, 'default': 35, 'help': 'Maximum peptide length. Default: 35.'},
    },
    'peptide_specificity': {
        'flags': ('--peptide-specificity',),
        'kwargs': {
            'default': 'full',
            'choices': ['full', 'n-specific', 'c-specific', 'semi'],
            'help': 'Cleavage specificity. Default: full.',
        },
    },
    'peptide_initiator_methionine': {
        'flags': ('--peptide-initiator-methionine',),
        'kwargs': {
            'default': 'both',
            'choices': ['both', 'remove', 'retain'],
            'help': 'N-terminal methionine handling. Default: both.',
        },
    },
    'peptide_isobaric': {
        'flags': ('--peptide-isobaric',),
        'kwargs': {
            'action': 'store_true',
            'help': 'Normalize I to L before peptide comparison and output.',
        },
    },
    'rebuild_peptide_sqlite': {
        'flags': ('--rebuild-peptide-sqlite',),
        'kwargs': {
            'action': 'store_true',
            'help': 'Atomically rebuild the peptide SQLite index before analysis.',
        },
    },
    'info_field': {
        'flags': ('--info_field',),
        'kwargs': {
            'help': 'fields to use in the INFO column of the vcf file to filter variants. Default None',
            'default': None,
        },
    },
    'info_field_thres': {
        'flags': ('--info_field_thres',),
        'kwargs': {
            'help': 'threhold for the info field. Default None, do not filter any variants. If set "--info_filed AF --info_field_thres 0.01", only keep variants with AF >= 0.01',
            'default': None,
        },
    },
    'chromosome': {
        'flags': ('-c', '--chromosome'),
        'kwargs': {
            'help': '''chromosome name/id, default="chr1" ''',
            'default': 'chr1',
            'type': str,
        },
    },
    'file_vcf': {
        'flags': ('-i', '--file_vcf'),
        'kwargs': {
            'help': '''input vcf file. It can be a gzip file. If multiple vcf files are provided, use "," to join the file names. For example, "--file_vcf file1.vcf,file2.vcf". A pattern match is also supported, but quote is required to get it work. For example '--file_vcf "file*.vcf"'. A TSV manifest with first header column "filepath" is also supported for one-VCF-per-sample population mode. Optional manifest columns: sample, name_use. ''',
            'required': True,
        },
    },
    'outprefix_vcf2mutation': {
        'flags': ('-o', '--outprefix'),
        'kwargs': {
            'help': 'output prefix to store the two output dataframes, default: None, do not write the result to files. file will be outprefix_1/2.tsv',
            'default': None,
        },
    },
    'protein_fasta': {
        'flags': ('-p', '--protein'),
        'kwargs': {
            'help': 'protein sequences in fasta format. It can be a gzip file.',
            'required': True,
        },
    },
    'outfile': {
        'flags': ('-o', '--outfile'),
        'kwargs': {
            'help': 'file to store the output PEFF file',
            'required': True,
        },
    },
    'mutation_annotation': {
        'flags': ('-m', '--mutation'),
        'kwargs': {
            'help': 'mutation annotation file generated by PrecisonProDB',
            'required': True,
        },
    },
    'peff_test': {
        'flags': ('-T', '--TEST'),
        'kwargs': {
            'help': 'whether to test the variant annotation is valid. Default: False',
            'action': 'store_true',
        },
    },
    'uniprots': {
        'flags': ('-u', '--uniprots'),
        'kwargs': {
            'help': 'uniprot sequences in fasta format. It can be a gzip file. If more than one files, join with ",". Default: None',
            'default': None,
        },
    },
    'changed': {
        'flags': ('-c', '--changed'),
        'kwargs': {
            'help': 'uniprot changed. A table generated by PrecisonProDB, including the relationship of uniprot_id and ensembl_id. It can be a gzip file.  Default: None',
            'default': None,
        },
    },
    'outfile2': {
        'flags': ('-O', '--outfile2'),
        'kwargs': {
            'help': 'file to store the output PEFF file based on uniprot sequences.  Default: None',
            'default': None,
        },
    },
    'files_uniprot': {
        'flags': ('-u', '--files_uniprot'),
        'kwargs': {
            'help': 'uniprot proteins. If more than one files, join by ","',
            'required': True,
        },
    },
    'files_ref': {
        'flags': ('-r', '--files_ref'),
        'kwargs': {
            'help': 'reference proteins to match with uniprot proteins. If more than one files, join by ","',
        },
    },
    'files_alt': {
        'flags': ('-a', '--files_alt'),
        'kwargs': {
            'help': 'altered reference proteins. If more than one files, join by ",". The order should be the same as files_ref',
        },
    },
    'outprefix_extract': {
        'flags': ('-o', '--outprefix'),
        'kwargs': {
            'help': 'prefix for output files. default:"perGeno"',
            'default': 'perGeno',
        },
    },
    'length_min': {
        'flags': ('-m', '--length_min'),
        'kwargs': {
            'help': 'minumum length required when matching UniProt sequences with sequences in files_ref. default: "20"',
            'default': 20,
            'type': int,
        },
    },
}


ARGUMENT_GROUPS = {
    'reference_inputs': ('genome', 'gtf', 'protein'),
    'variant_input_detailed': ('mutations_detailed',),
    'variant_input_basic': ('mutations_basic',),
    'variant_input_perchrom': ('mutations_perchrom',),
    'annotation': ('datatype', 'protein_keyword'),
    'datatype_only': ('datatype',),
    'runtime_output': ('threads', 'out'),
    'vcf_selection': ('no_filter', 'sample_vcf', 'all_chromosomes'),
    'vcf_selection_mixed': ('no_filter', 'sample_mixed', 'all_chromosomes'),
    'vcf_info_filters': ('info_field', 'info_field_thres'),
    'download_uniprot': ('download', 'uniprot', 'uniprot_min_len', 'peff'),
    'cleanup': ('keep_all',),
    'sqlite': ('sqlite',),
    'peptide': (
        'peptide', 'peptide_sqlite', 'peptide_enzyme', 'peptide_missed_cleavages',
        'peptide_min_length', 'peptide_max_length', 'peptide_specificity',
        'peptide_initiator_methionine', 'peptide_isobaric', 'rebuild_peptide_sqlite',
    ),
    'perchrom': ('chromosome',),
    'perchrom_sqlite': ('chromosome', 'sample_perchrom'),
    'vcf_conversion': (
        'file_vcf',
        'outprefix_vcf2mutation',
        'sample_vcf',
        'no_filter',
        'all_chromosomes',
        'info_field',
        'info_field_thres',
        'threads_vcf2mutation',
    ),
    'peff_generation': ('protein_fasta', 'outfile', 'mutation_annotation', 'peff_test', 'uniprots', 'changed', 'outfile2'),
    'uniprot_matching': ('files_uniprot', 'files_ref', 'files_alt', 'outprefix_extract', 'length_min'),
}


def add_argument_set(parser, *group_names, overrides=None):
    overrides = overrides or {}
    added = set()
    for group_name in group_names:
        for spec_name in ARGUMENT_GROUPS[group_name]:
            if spec_name in added:
                continue
            spec = ARGUMENT_SPECS[spec_name]
            override = dict(overrides.get(spec_name, {}))
            flags = tuple(override.pop('flags', spec['flags']))
            kwargs = dict(spec['kwargs'])
            kwargs.update(override)
            parser.add_argument(*flags, **kwargs)
            added.add(spec_name)
    return parser
