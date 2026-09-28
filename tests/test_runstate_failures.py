"""Regression tests for safe reuse and complete reference-only output."""

import sqlite3
import subprocess
import sys

import pandas as pd
import pytest

from precisionprodb import PrecisionProDB_Sqlite as sqlite_pipeline
from precisionprodb import buildSqlite
from precisionprodb.PrecisionProDB_core import PerGeno
from precisionprodb.runstate import RunState, sqlite_owned_by_run
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
    haplotype_temp = tmp_path / 'sample_1_temp'
    haplotype_temp.mkdir()
    marker = haplotype_temp / '1.perChromFinished'
    marker.write_text('1')
    first.complete()
    assert state().prepare() is True
    marker.write_text('stale')
    with pytest.raises(ValueError, match='intermediate files changed'):
        state().prepare()
    assert state(force=True).prepare() is False
    assert not haplotype_temp.exists()


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
    with pytest.raises(ValueError, match='inputs/settings changed'):
        convertVCF2MutationComplex(str(source), prefix, individual_input='B')
    assert convertVCF2MutationComplex(str(source), prefix, individual_input='B', force=True) == ['B__1', 'B__2']
    assert list((tmp_path / 'converted.tsv.archive').glob('*/converted.tsv'))
