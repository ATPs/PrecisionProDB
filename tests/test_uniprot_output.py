import csv

from Bio import SeqIO

from precisionprodb.extractMutatedUniprot import extractMutatedUniprot
from precisionprodb.PrecisionProDB_test import get_cmd_to_run, validate_case_outputs


def write_fasta(path, entries):
    with open(path, 'w') as handle:
        for record_id, sequence, description in entries:
            suffix = f' {description}' if description else ''
            handle.write(f'>{record_id}{suffix}\n{sequence}\n')


def read_fasta(path):
    return [(record.id, str(record.seq), record.description) for record in SeqIO.parse(path, 'fasta')]


def test_duplicate_reference_ids_and_alternatives_are_preserved_stably(tmp_path):
    uniprot = tmp_path / 'uniprot.fa'
    ref = tmp_path / 'ref.fa'
    ref_reversed = tmp_path / 'ref_reversed.fa'
    alt = tmp_path / 'alt.fa'
    write_fasta(uniprot, [('UP1', 'MABCDE', 'UniProt description')])
    reference_entries = [
        ('R2', 'MABCDE', ''),
        ('R1', 'MABCDE', ''),
    ]
    write_fasta(ref, reference_entries)
    write_fasta(ref_reversed, list(reversed(reference_entries)))
    write_fasta(alt, [
        ('R2__1', 'MABYDE', 'individual2\tchanged'),
        ('R1__1', 'MABXDE', 'individual1\tchanged'),
        ('R2__2', 'MABXDE', 'individual3\tchanged'),
        ('R2__3', 'MABCDE', 'reference\tunchanged'),
    ])

    first_prefix = str(tmp_path / 'first')
    reversed_prefix = str(tmp_path / 'reversed')
    extractMutatedUniprot(str(uniprot), str(ref), str(alt), first_prefix, length_min=1)
    extractMutatedUniprot(
        str(uniprot), str(ref_reversed), str(alt), reversed_prefix, length_min=1
    )

    for extension in ('.uniprot_changed.tsv', '.uniprot_changed.fa', '.uniprot_all.fa'):
        assert (tmp_path / f'first{extension}').read_bytes() == (
            tmp_path / f'reversed{extension}'
        ).read_bytes()

    with open(first_prefix + '.uniprot_changed.tsv', newline='') as handle:
        relationships = list(csv.DictReader(handle, delimiter='\t'))
    assert relationships == [
        {'uniprot_id': 'UP1', 'ref_id': 'R1'},
        {'uniprot_id': 'UP1', 'ref_id': 'R2'},
    ]

    changed = read_fasta(first_prefix + '.uniprot_changed.fa')
    assert [(record_id, sequence) for record_id, sequence, _ in changed] == [
        ('UP1__1', 'MABXDE'),
        ('UP1__2', 'MABYDE'),
    ]
    all_records = read_fasta(first_prefix + '.uniprot_all.fa')
    assert [(record_id, sequence) for record_id, sequence, _ in all_records] == [
        ('UP1__reference', 'MABCDE'),
        ('UP1__1', 'MABXDE'),
        ('UP1__2', 'MABYDE'),
    ]


def test_unmatched_or_unaltered_uniprot_remains_in_all_fasta(tmp_path):
    uniprot = tmp_path / 'uniprot.fa'
    ref = tmp_path / 'ref.fa'
    alt = tmp_path / 'alt.fa'
    write_fasta(uniprot, [
        ('UP1', 'MABCDE', 'matched'),
        ('UP2', 'MQQQQQ', 'unmatched'),
    ])
    write_fasta(ref, [('R1', 'MABCDE', '')])
    write_fasta(alt, [('R1', 'MABCDE', 'reference\tunchanged')])

    prefix = str(tmp_path / 'result')
    extractMutatedUniprot(str(uniprot), str(ref), str(alt), prefix, length_min=1)

    assert read_fasta(prefix + '.uniprot_changed.fa') == []
    assert read_fasta(prefix + '.uniprot_all.fa') == [
        ('UP1', 'MABCDE', 'UP1 matched\tunchanged'),
        ('UP2', 'MQQQQQ', 'UP2 unmatched\tunchanged'),
    ]


def test_regression_driver_disables_sqlite_with_explicit_none(tmp_path):
    root = tmp_path / 'repo'
    output = tmp_path / 'output'
    commands = get_cmd_to_run(
        'Ensembl',
        'tsv',
        'no_sqlite',
        str(output),
        {'tsv': '/input/variants.tsv'},
        {'Ensembl': {
            'genome': '/input/genome.fa',
            'protein': '/input/protein.fa',
            'gtf': '/input/genes.gtf',
            'datatype': 'Ensembl_GTF',
        }},
        str(root),
    )
    assert len(commands) == 1
    command = commands[0]
    assert command[0]
    assert command[command.index('--sqlite') + 1] == 'NONE'
    assert command[command.index('-m') + 1] == '/input/variants.tsv'


def test_regression_driver_validates_mutation_and_fasta_identifiers(tmp_path):
    prefix = str(tmp_path / 'run')
    (tmp_path / 'run.pergeno.aa_mutations.csv').write_text(
        'protein_id\tprotein_id_fasta\nP1\tP1__1\n'
    )
    write_fasta(tmp_path / 'run.pergeno.protein_changed.fa', [('P1__1', 'MABCDE', 'changed')])
    write_fasta(tmp_path / 'run.pergeno.protein_all.fa', [('P1__1', 'MABCDE', 'changed')])
    write_fasta(tmp_path / 'run.pergeno.protein_PEFF.fa', [('P1__1', 'MABCDE', 'PEFF')])

    validate_case_outputs(prefix, 'Ensembl')

    write_fasta(tmp_path / 'run.pergeno.protein_all.fa', [])
    try:
        validate_case_outputs(prefix, 'Ensembl')
    except RuntimeError as error:
        assert 'absent from all-protein FASTA' in str(error)
    else:
        raise AssertionError('expected changed/all FASTA mismatch to fail validation')
