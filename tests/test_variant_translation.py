import sqlite3

import pandas as pd
import pytest

from precisionprodb import buildSqlite, perChrom, perChromSqlite


@pytest.fixture(autouse=True)
def _reset_translation_caches():
    perChrom.reset_translation_caches()
    yield
    perChrom.reset_translation_caches()


def _translation_row(sequence, strand='+'):
    return pd.Series(
        {
            'seqname': '1',
            'genomicLocs': [(0, len(sequence))],
            'CDSplus': sequence,
            'AA_seq': 'MTK',
            'AA_translate': 'MTK',
            'frame': 0,
            'strand': strand,
            'mutations': [0],
        },
        name='protein1',
    )


def test_sqlite_interval_lookup_uses_1_based_inclusive_boundaries():
    con = sqlite3.connect(':memory:')
    con.execute(
        'CREATE TABLE genomicLocs_1 '
        '(protein_id TEXT, genomicLocs_start INTEGER, genomicLocs_end INTEGER)'
    )
    con.execute('INSERT INTO genomicLocs_1 VALUES (?, ?, ?)', ('protein1', 101, 110))
    con.commit()

    assert buildSqlite.get_protein_id_from_genomicLocs(con, '1', 101) == ['protein1']
    assert buildSqlite.get_protein_id_from_genomicLocs(con, '1', 110) == ['protein1']
    assert buildSqlite.get_protein_id_from_genomicLocs(con, '1', 100) == []
    assert buildSqlite.get_protein_id_from_genomicLocs(con, '1', 111) == []
    assert buildSqlite.get_protein_id_from_genomicLocs(con, '1', 101, 110) == ['protein1']
    assert buildSqlite.get_protein_id_from_genomicLocs(con, '1', 101, 111) == []

    queries = [('1', 101, 101), ('1', 110, 110), ('1', 100, 100), ('1', 101, 110)]
    assert buildSqlite.get_protein_id_from_genomicLocs(
        con, [query[0] for query in queries], [query[1] for query in queries],
        [query[2] for query in queries], threads=1,
    ) == [['protein1'], ['protein1'], [], ['protein1']]
    con.close()


def test_dataframe_interval_helpers_match_scalar_and_batch_boundaries():
    intervals = pd.DataFrame(
        {
            'protein_id': ['protein1'],
            'genomicLocs_start': [101],
            'genomicLocs_end': [110],
        }
    )

    assert buildSqlite.get_protein_ids_from_tdf(intervals, 101) == ['protein1']
    assert buildSqlite.get_protein_ids_from_tdf(intervals, 100) == []
    assert buildSqlite.get_protein_ids_from_tdf_bisect(intervals, 110) == ['protein1']
    assert buildSqlite.get_protein_ids_from_tdf_batch(
        intervals, [(101, 101), (110, 110), (100, 100), (101, 111)]
    ) == [['protein1'], ['protein1'], [], []]


def test_multibase_substitution_translates_all_alt_bases_from_string():
    mutations = perChrom.parse_mutation('1-4-AC-GT', chromosome='1')
    row = _translation_row('ATGACCAAATAA')

    result = perChrom.translateCDSplusWithMut2(row, mutations)

    assert result['new_AA'] == 'MVK'
    assert result['variant_AA'] == 'T2V(1-4-AC-GT)'


def test_multibase_substitution_is_oriented_on_reverse_strand():
    # The forward reference GT at positions 8-9 reverse-complements to AC in CDS.
    mutations = perChrom.parse_mutation('1-8-GT-AC', chromosome='1')
    row = _translation_row('TTATTTGGTCAT', strand='-')

    result = perChrom.translateCDSplusWithMut2(row, mutations)

    assert result['new_AA'] == 'MVK'


def test_full_reference_allele_is_checked_before_applying_mnv():
    mutations = perChrom.parse_mutation('1-4-TC-GT', chromosome='1')
    row = _translation_row('ATGACCAAATAA')

    with pytest.raises(RuntimeError, match='reference mismatch'):
        perChrom.translateCDSplusWithMut2(row, mutations)


@pytest.mark.parametrize(
    ('variant', 'expected_cds'),
    [
        ('1-4-AC-A', 'ATGA CAAATAA'.replace(' ', '')),
        ('1-4-A-AT', 'ATGATCCAAATAA'),
        ('1-4-AC-GTT', 'ATGGTTCAAATAA'),
    ],
)
def test_length_changing_alleles_apply_as_complete_transcript_edits(variant, expected_cds):
    mutations = perChrom.parse_mutation(variant, chromosome='1')
    row = _translation_row('ATGACCAAATAA')
    cds, _, _, _, _ = perChrom.create_df_CDSplus_for_transcript_id(row)
    cds.attrs['transcript_id'] = row.name
    edits = perChrom.getMut_helper([0], '+', mutations, cds)
    merged = perChrom._merge_CDSplus_with_mutations(row.name, cds, edits)
    altered = merged['alt'].where(merged['alt'].notna(), merged['bases'])

    assert ''.join(altered) == expected_cds


def test_empty_translation_outputs_keep_headers_and_empty_fasta(tmp_path):
    outprefix = str(tmp_path / 'empty')

    perChrom.save_mutation_and_proteins(pd.DataFrame(), outprefix)

    table = pd.read_csv(outprefix + '.aa_mutations.csv', sep='\t')
    assert table.empty
    assert 'protein_id' in table.columns
    assert (tmp_path / 'empty.mutated_protein.fa').read_text() == ''


@pytest.mark.parametrize('protein_id', ['P', 'P__1'])
def test_standard_output_ids_match_fasta_ids_exactly(tmp_path, protein_id):
    outprefix = str(tmp_path / 'changed')
    results = pd.DataFrame(
        [{
            'protein_id_fasta': protein_id,
            'protein_description': protein_id + ' description',
            'AA_seq': 'MAK',
            'new_AA': 'MVK',
        }],
        index=pd.Index(['canonical'], name='protein_id'),
    )

    perChrom.save_mutation_and_proteins(results, outprefix)

    table = pd.read_csv(outprefix + '.aa_mutations.csv', sep='\t')
    assert table['protein_id_fasta'].tolist() == [protein_id]
    fasta_header = (tmp_path / 'changed.mutated_protein.fa').read_text().splitlines()[0]
    assert fasta_header.split()[0] == '>' + protein_id


def test_sqlite_save_handles_missing_new_sequence_column(tmp_path):
    outprefix = str(tmp_path / 'empty_sqlite')

    perChromSqlite.save_mutation_and_proteins(
        pd.DataFrame(columns=['individual']), outprefix
    )

    table = pd.read_csv(outprefix + '.aa_mutations.csv', sep='\t')
    assert table.empty
    assert 'individual' in table.columns
    assert (tmp_path / 'empty_sqlite.mutated_protein.fa').read_text() == ''
