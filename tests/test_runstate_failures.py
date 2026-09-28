"""Regression tests for safe reuse and complete reference-only output."""

import gzip
import json
import sqlite3
import subprocess
import sys

import pandas as pd
import pytest

from precisionprodb import PrecisionProDB_Sqlite as sqlite_pipeline
from precisionprodb import buildSqlite
from precisionprodb.PrecisionProDB_core import PerGeno
from precisionprodb.runstate import RunState, mutation_files, sqlite_owned_by_run
from precisionprodb.vcf2mutation import convertVCF2MutationComplex


def test_runstate_refuses_changed_input_and_force_rebuilds(tmp_path):
    variants = tmp_path / 'variants.tsv'
    variants.write_text('chr\tpos\tref\talt\n1\t10\tA\tG\n')
    prefix = str(tmp_path / 'sample')

    def state(force=False):
        return RunState.from_inputs(
            prefix, str(variants), '', '', '', '', '', False,
            {'sample': 'A'}, force=force, keep_all=True,
        )

    first = state()
    first.prepare()
    (tmp_path / 'sample.pergeno.protein_all.fa').write_text('>P\nMAK\n')
    (tmp_path / 'sample.pergeno.protein_changed.fa').write_text('')
    (tmp_path / 'sample.pergeno.aa_mutations.csv').write_text('protein_id_fasta\n')
    first.complete()
    assert state().prepare() is True
    variants.write_text('chr\tpos\tref\talt\n1\t20\tC\tT\n')
    with pytest.raises(ValueError, match='inputs/settings changed'):
        state().prepare()
    rebuilt = state(force=True)
    rebuilt.prepare()
    assert not (tmp_path / 'sample.pergeno.protein_all.fa').exists()
    assert list((tmp_path / 'sample.archive').glob('*/sample.pergeno.protein_all.fa'))


def test_generated_explicit_sqlite_is_stable_on_same_prefix_reuse(tmp_path):
    variants = tmp_path / 'variants.tsv'
    genome = tmp_path / 'genome.fa'
    gtf = tmp_path / 'genes.gtf'
    protein = tmp_path / 'protein.fa'
    for path in (variants, genome, gtf, protein):
        path.write_text('fixture\n')
    prefix = str(tmp_path / 'sample')
    database = str(tmp_path / 'annotation.sqlite')

    def state(force=False):
        owned = sqlite_owned_by_run(prefix, database, str(genome), str(gtf), str(protein))
        return RunState.from_inputs(
            prefix, str(variants), database, str(genome), str(gtf),
            str(protein), '', False,
            {'sqlite': str(database), 'owned_sqlite': owned},
            force=force, owned_sqlite=owned,
        )

    first = state()
    assert first.owned_sqlite
    assert first.prepare() is False
    (tmp_path / 'annotation.sqlite').write_text('built annotation')
    (tmp_path / 'sample.pergeno.protein_all.fa').write_text('>P\nMAK\n')
    (tmp_path / 'sample.pergeno.protein_changed.fa').write_text('')
    (tmp_path / 'sample.pergeno.aa_mutations.csv').write_text('protein_id_fasta\n')
    first.complete()
    assert state().prepare() is True
    (tmp_path / 'annotation.sqlite').write_text('modified annotation')
    with pytest.raises(ValueError, match='generated output changed'):
        state().prepare()


def test_prebuilt_explicit_sqlite_is_external_even_at_default_path(tmp_path):
    prefix = str(tmp_path / 'sample')
    database = prefix + '.sqlite'
    (tmp_path / 'sample.sqlite').write_text('prebuilt annotation')
    assert not sqlite_owned_by_run(
        prefix, database, 'genome.fa', 'genes.gtf', 'protein.fa',
        default_mode=False,
    )


@pytest.mark.parametrize('keep_all', [False, True])
@pytest.mark.parametrize(
    ('receipt_version', 'recorded_owned'),
    [(2, False), (2, True), (3, False), (3, True), (None, None)],
)
def test_force_preserves_unverified_default_sqlite(
        tmp_path, keep_all, receipt_version, recorded_owned):
    prefix = str(tmp_path / 'sample')
    database = tmp_path / 'sample.sqlite'
    original = b'annotation database to preserve'
    database.write_bytes(original)
    if receipt_version is not None:
        (tmp_path / 'sample.run.json').write_text(json.dumps({
            'format_version': receipt_version,
            'status': 'complete',
            'settings': {'sqlite': str(database), 'owned_sqlite': recorded_owned},
        }))
    output = tmp_path / 'sample.pergeno.protein_all.fa'
    output.write_text('>P\nMAK\n')

    owned = sqlite_owned_by_run(
        prefix, str(database), 'genome.fa', 'genes.gtf', 'protein.fa',
        default_mode=True,
    )
    expected_owned = receipt_version == 3 and recorded_owned is True
    assert owned is expected_owned
    state = RunState(
        prefix, [], {'sqlite': str(database), 'owned_sqlite': owned},
        force=True, keep_all=keep_all, owned_sqlite=owned,
    )
    assert state.prepare() is False
    assert not output.exists()
    if expected_owned:
        assert not database.exists()
        if keep_all:
            archived = list((tmp_path / 'sample.archive').glob('*/sample.sqlite'))
            assert len(archived) == 1
            assert archived[0].read_bytes() == original
    else:
        assert database.read_bytes() == original
        assert not list((tmp_path / 'sample.archive').glob('*/sample.sqlite'))


def test_external_sqlite_still_fingerprints_explicit_protein_fasta(tmp_path):
    variants = tmp_path / 'variants.tsv'
    variants.write_text('chr\tpos\tref\talt\n')
    database = tmp_path / 'annotation.sqlite'
    database.write_bytes(b'external annotation')
    protein = tmp_path / 'proteins.fa'
    protein.write_text('>P\nMAK\n')
    prefix = str(tmp_path / 'sample')

    def state():
        return RunState.from_inputs(
            prefix, str(variants), str(database), '', '', str(protein), '',
            False, {'sqlite': str(database)},
        )

    first = state()
    first.prepare()
    (tmp_path / 'sample.pergeno.protein_all.fa').write_text('>P\nMAK\n')
    (tmp_path / 'sample.pergeno.protein_changed.fa').write_text('')
    (tmp_path / 'sample.pergeno.aa_mutations.csv').write_text('protein_id_fasta\n')
    first.complete()
    protein.write_text('>P\nMVK\n')
    with pytest.raises(ValueError, match='inputs/settings changed'):
        state().prepare()


def test_retained_legacy_haplotype_intermediates_are_checked(tmp_path):
    variants = tmp_path / 'variants.vcf'
    variants.write_text('fixture\n')
    prefix = str(tmp_path / 'sample')

    def state(force=False):
        return RunState.from_inputs(
            prefix, str(variants), '', '', '', '', '', False,
            {'sample': 'A'}, force=force, keep_all=True,
        )

    first = state()
    assert first.prepare() is False
    (tmp_path / 'sample.pergeno.protein_all.fa').write_text('>P\nMAK\n')
    (tmp_path / 'sample.pergeno.protein_changed.fa').write_text('')
    (tmp_path / 'sample.pergeno.aa_mutations.csv').write_text('protein_id_fasta\n')
    independent_fasta = tmp_path / 'sample_1.pergeno.protein_all.fa'
    independent_fasta.write_text('>independent\nMAK\n')
    haplotype_temp = tmp_path / 'sample_temp' / 'haplotypes'
    haplotype_temp.mkdir(parents=True)
    marker = haplotype_temp / '1.perChromFinished'
    marker.write_text('1')
    first.complete()
    assert state().prepare() is True
    marker.write_text('stale')
    with pytest.raises(ValueError, match='intermediate files changed'):
        state().prepare()
    assert state(force=True).prepare() is False
    assert not haplotype_temp.exists()
    assert independent_fasta.read_text() == '>independent\nMAK\n'


def test_version_two_receipt_requires_explicit_rebuild(tmp_path):
    prefix = str(tmp_path / 'sample')
    database = tmp_path / 'external.sqlite'
    database.write_bytes(b'external annotation')
    (tmp_path / 'sample.run.json').write_text(json.dumps({
        'format_version': 2,
        'status': 'complete',
        'settings': {'sqlite': str(database), 'owned_sqlite': True},
    }))
    (tmp_path / 'sample.pergeno.protein_all.fa').write_text('>P\nMAK\n')
    neighbor = tmp_path / 'sample_1.pergeno.protein_all.fa'
    neighbor.write_text('>neighbor\nMST\n')

    state = RunState(
        prefix, [], {'sqlite': str(database)}, owned_sqlite=False,
    )
    with pytest.raises(ValueError, match='receipt format 2'):
        state.prepare()
    rebuilt = RunState(
        prefix, [], {'sqlite': str(database)}, force=True,
        owned_sqlite=False,
    )
    assert rebuilt.prepare() is False
    assert database.read_bytes() == b'external annotation'
    assert neighbor.read_text() == '>neighbor\nMST\n'


def test_gzipped_manifest_inputs_are_fingerprinted(tmp_path):
    vcf = tmp_path / 'sample.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n')
    manifest = tmp_path / 'manifest.tsv.gz'
    with gzip.open(manifest, 'wt') as handle:
        handle.write(f'filepath\n  {vcf}  \n')

    assert mutation_files(str(manifest), is_manifest=True) == [
        str(manifest), str(vcf)
    ]
    state = RunState.from_inputs(
        str(tmp_path / 'out'), str(manifest), '', '', '', '', '', True, {},
    )
    assert [item['path'] for item in state.inputs] == [
        str(manifest.resolve()), str(vcf.resolve())
    ]


def test_receipt_writer_creates_nested_output_directory(tmp_path):
    prefix = str(tmp_path / 'new' / 'nested' / 'sample')
    state = RunState(prefix, [], {})
    assert state.prepare() is False
    assert (tmp_path / 'new' / 'nested' / 'sample.run.json').is_file()


def test_standalone_sqlite_build_rejects_stale_split_files(tmp_path):
    prefix = str(tmp_path / 'build')
    intermediates = tmp_path / 'build_temp'
    intermediates.mkdir()
    (intermediates / 'split.done').write_text('old')
    with pytest.raises(ValueError, match='unverified SQLite build intermediates'):
        buildSqlite.create_sqlite(
            str(tmp_path / 'annotation.sqlite'), '', '', '', prefix,
            'gtf', 'auto',
        )


def test_empty_sqlite_mutations_produce_reference_only_fasta(tmp_path, monkeypatch):
    source = tmp_path / 'ref.sqlite'
    with sqlite3.connect(source) as connection:
        connection.execute(
            'CREATE TABLE protein_description '
            '(protein_id_fasta TEXT, protein_description TEXT, AA_seq TEXT)'
        )
        connection.execute("INSERT INTO protein_description VALUES ('P', 'P', 'MAK')")
    variants = tmp_path / 'empty.tsv'
    variants.write_text('chr\tpos\tref\talt\n')
    monkeypatch.setattr(PerGeno, 'splitMutationByChromosomeLarge', lambda *a, **k: [])
    prefix = str(tmp_path / 'out')
    sqlite_pipeline.runPerChomSqlite(
        str(source), str(variants), 1, prefix, 'auto', 'gtf', True,
        '', ['1'], ['1 chromosome 1,'], '',
    )
    assert (tmp_path / 'out.pergeno.protein_all.fa').read_text() == '>P\tunchanged\nMAK\n'
    assert (tmp_path / 'out.pergeno.protein_changed.fa').read_text() == ''
    table = pd.read_csv(tmp_path / 'out.pergeno.aa_mutations.csv', sep='\t')
    assert table.empty
    assert 'protein_id' in table.columns


def test_direct_variant_database_includes_changed_and_reference_proteins(tmp_path):
    source = tmp_path / 'ref.sqlite'
    with sqlite3.connect(source) as connection:
        connection.execute(
            'CREATE TABLE protein_description (protein_description TEXT, AA_seq TEXT)'
        )
        connection.executemany(
            'INSERT INTO protein_description VALUES (?, ?)',
            [('P1 description', 'MAK'), ('P2 description', 'MST')],
        )
    prefix = str(tmp_path / 'out')
    (tmp_path / 'out.pergeno.mutated_protein.fa').write_text(
        '>P1 description\tchanged\nMVK\n'
    )
    sqlite_pipeline.write_direct_variant_protein_database(str(source), prefix)
    assert (tmp_path / 'out.pergeno.protein_changed.fa').read_text() == (
        '>P1 description\tchanged\nMVK\n'
    )
    all_fasta = (tmp_path / 'out.pergeno.protein_all.fa').read_text()
    assert '>P1 description\tchanged\nMVK\n' in all_fasta
    assert '>P2 description\tunchanged\nMST\n' in all_fasta
    assert '>P1 description\tunchanged' not in all_fasta


def test_chromosome_failure_propagates(tmp_path, monkeypatch):
    class FailingChromosome:
        def __init__(self, **kwargs):
            pass

        def run_perChrom(self):
            raise RuntimeError('simulated failure')

    monkeypatch.setattr(sqlite_pipeline.perChromSqlite, 'PerChrom_sqlite', FailingChromosome)
    with pytest.raises(RuntimeError, match='failed processing chromosome chr1'):
        sqlite_pipeline.runSinglePerChromSqlite(
            'unused', 'unused', str(tmp_path), 1, 'chr1', 'gtf', ''
        )


def test_invalid_sqlite_exits_nonzero(tmp_path):
    source = tmp_path / 'bad.sqlite'
    source.write_text('not a database')
    variants = tmp_path / 'mutations.tsv'
    variants.write_text('chr\tpos\tref\talt\n')
    result = subprocess.run([
        sys.executable, '-m', 'precisionprodb.PrecisionProDB',
        '-S', str(source), '-m', str(variants), '-o', str(tmp_path / 'out'),
    ], capture_output=True, text=True)
    assert result.returncode == 2
    assert 'invalid annotation SQLite' in result.stderr


def test_converter_cache_checks_selection_and_force(tmp_path):
    source = tmp_path / 'two_samples.vcf'
    source.write_text(
        '##fileformat=VCFv4.2\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n'
        '1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1\n'
    )
    prefix = str(tmp_path / 'converted')
    assert convertVCF2MutationComplex(str(source), prefix, individual_input='A') == ['A__1', 'A__2']
    receipt = json.loads((tmp_path / 'converted.tsv.cache.json').read_text())
    assert receipt['format_version'] == 2
    with pytest.raises(ValueError, match='inputs/settings changed'):
        convertVCF2MutationComplex(str(source), prefix, individual_input='B')
    assert convertVCF2MutationComplex(str(source), prefix, individual_input='B', force=True) == ['B__1', 'B__2']
    assert list((tmp_path / 'converted.tsv.archive').glob('*/converted.tsv'))
