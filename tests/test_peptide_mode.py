import csv
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import subprocess
import sys

import pandas as pd
import pytest

from precisionprodb.peptide import (
    iter_alt_specific_peptides,
    merge_novel_peptide_parts,
    write_novel_peptide_part,
)
from precisionprodb.peptideSqlite import (
    KnownPeptideIndex,
    PeptideConfig,
    build_peptide_sqlite,
    default_peptide_sqlite_path,
    ensure_peptide_sqlite,
    validate_peptide_sqlite,
)


def make_annotation_sqlite(path):
    connection = sqlite3.connect(path)
    connection.execute(
        'CREATE TABLE protein_description (protein_id TEXT PRIMARY KEY, AA_seq TEXT)'
    )
    connection.executemany(
        'INSERT INTO protein_description VALUES (?, ?)',
        [
            ('P1', 'MPEPTIDEKAAAAK'),
            ('P2', 'AAATK'),
            ('P3', 'AKRPQK'),
        ],
    )
    connection.commit()
    connection.close()


def update_metadata(path, key, value):
    connection = sqlite3.connect(path)
    connection.execute('UPDATE metadata SET value = ? WHERE key = ?', (value, key))
    connection.commit()
    connection.close()


@pytest.fixture()
def index(tmp_path):
    source = tmp_path / 'annotation.sqlite'
    output = tmp_path / 'known.peptides.sqlite'
    make_annotation_sqlite(source)
    config = PeptideConfig(
        min_length=3,
        max_length=30,
        missed_cleavages=0,
        initiator_methionine='retain',
    )
    metadata = build_peptide_sqlite(source, output, config)
    return source, output, config, metadata


def test_index_build_lookup_and_validation(index):
    source, output, config, metadata = index
    assert metadata['build_complete'] == '1'
    assert int(metadata['known_peptide_count']) > 0
    validated = validate_peptide_sqlite(output, config, source)
    assert validated['digestion_config_hash'] == config.config_hash
    with KnownPeptideIndex(output) as known:
        assert known.lookup_many(['MPEPTIDEK', 'AAATK', 'QQQK']) == {
            'MPEPTIDEK',
            'AAATK',
        }
        assert known.lookup_many(['QQQK']) == set()


def test_configuration_mismatch_is_rejected(index):
    source, output, _config, _metadata = index
    with pytest.raises(ValueError, match='configuration'):
        validate_peptide_sqlite(output, PeptideConfig(min_length=4, max_length=30), source)


def test_existing_mismatch_requires_explicit_rebuild(index):
    source, output, _config, _metadata = index
    new_config = PeptideConfig(
        min_length=4,
        max_length=30,
        missed_cleavages=0,
        initiator_methionine='retain',
    )
    with pytest.raises(ValueError, match='configuration'):
        ensure_peptide_sqlite(source, new_config, output)
    assert ensure_peptide_sqlite(source, new_config, output, rebuild=True) == str(output)
    validated = validate_peptide_sqlite(output, new_config, source)
    assert validated['digestion_config_hash'] == new_config.config_hash


def test_selector_and_rebuild_require_explicit_peptide_gate(tmp_path):
    selector = str(tmp_path / 'known.peptides.sqlite')
    selector_result = subprocess.run(
        [
            sys.executable,
            '-m',
            'precisionprodb.PrecisionProDB',
            '-m',
            'ignored',
            '--peptide-sqlite',
            selector,
        ],
        capture_output=True,
        text=True,
    )
    rebuild_result = subprocess.run(
        [
            sys.executable,
            '-m',
            'precisionprodb.PrecisionProDB',
            '-m',
            'ignored',
            '--rebuild-peptide-sqlite',
        ],
        capture_output=True,
        text=True,
    )
    assert selector_result.returncode == 2
    assert rebuild_result.returncode == 2
    assert 'require --peptide' in selector_result.stderr
    assert 'require --peptide' in rebuild_result.stderr


def test_novel_writer_filters_reference_and_global_known(index, tmp_path):
    _source, output, config, _metadata = index
    changed = pd.DataFrame([
        {
            'protein_id': 'P1',
            'protein_id_fasta': 'P1__1',
            'AA_seq': 'MPEPTIDEKAAAAK',
            'new_AA': 'MPEPTIDEKAAATK',
            'individual': 'sampleA',
            'seqname': 'chr1',
            'variant_AA': 'A13T',
            'insertion_AA': '',
            'deletion_AA': '',
            'frameChange': False,
            'stopGain': False,
            'AA_stopGain': '',
            'stopLoss': False,
            'stopLoss_pos': '',
        },
        {
            'protein_id': 'P1',
            'protein_id_fasta': 'P1__2',
            'AA_seq': 'MPEPTIDEKAAAAK',
            'new_AA': 'MPEPTIDEKQQQK',
            'individual': 'sampleB',
            'seqname': 'chr1',
            'variant_AA': 'AAAA>QQQ',
            'insertion_AA': '',
            'deletion_AA': '',
            'frameChange': False,
            'stopGain': False,
            'AA_stopGain': '',
            'stopLoss': False,
            'stopLoss_pos': '',
        },
    ])
    part = tmp_path / 'part.tsv'
    with KnownPeptideIndex(output) as known:
        write_novel_peptide_part(changed, known, part, config)
    with open(part, newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert [row['peptide_sequence'] for row in rows] == ['QQQK']
    assert rows[0]['protein_id_fasta'] == 'P1__2'
    final_tsv, final_fasta = merge_novel_peptide_parts([part], str(tmp_path / 'result'))
    assert 'PPDBpep_000001' in Path(final_tsv).read_text()
    assert Path(final_fasta).read_text() == '>PPDBpep_000001|n_mappings=1\nQQQK\n'


def test_parallel_extraction_matches_serial_output(index, tmp_path):
    _source, output, config, _metadata = index
    changed = pd.DataFrame([
        {
            'protein_id': 'P1',
            'protein_id_fasta': 'P1__1',
            'AA_seq': 'MPEPTIDEKAAAAK',
            'new_AA': 'MPEPTIDEKQQQK',
            'individual': 'sampleA',
            'seqname': 'chr1',
            'variant_AA': 'AAAA>QQQ',
            'insertion_AA': '',
            'deletion_AA': '',
            'frameChange': False,
            'stopGain': False,
            'AA_stopGain': '',
            'stopLoss': False,
            'stopLoss_pos': '',
        },
        {
            'protein_id': 'P3',
            'protein_id_fasta': 'P3__1',
            'AA_seq': 'AKRPQK',
            'new_AA': 'AKRPQW',
            'individual': 'sampleB',
            'seqname': 'chr2',
            'variant_AA': 'K6W',
            'insertion_AA': '',
            'deletion_AA': '',
            'frameChange': False,
            'stopGain': False,
            'AA_stopGain': '',
            'stopLoss': False,
            'stopLoss_pos': '',
        },
    ])
    serial = tmp_path / 'serial.tsv'
    parallel = tmp_path / 'parallel.tsv'
    with KnownPeptideIndex(output) as known:
        write_novel_peptide_part(changed, known, serial, config, threads=1)
    with KnownPeptideIndex(output) as known:
        write_novel_peptide_part(changed, known, parallel, config, threads=2)
    assert serial.read_bytes() == parallel.read_bytes()


def test_failed_build_does_not_publish_index(tmp_path):
    output = tmp_path / 'known.peptides.sqlite'
    with pytest.raises(FileNotFoundError):
        build_peptide_sqlite(tmp_path / 'missing.sqlite', output, PeptideConfig())
    assert not output.exists()
    assert not list(tmp_path.glob('known.peptides.sqlite.building.*'))


def test_sequence_subtraction_handles_cleavage_indel_and_stop_changes():
    config = PeptideConfig(
        min_length=3,
        max_length=30,
        missed_cleavages=0,
        initiator_methionine='retain',
    )
    cleavage_gain = list(iter_alt_specific_peptides('AAAKAAAAK', 'AAARAAAAK', config))
    deletion_junction = list(iter_alt_specific_peptides('AAAKVVVVK', 'AAAKVVK', config))
    stop_gain = list(iter_alt_specific_peptides('MPEPTIDEKQQQK', 'MPEPTIDEKQQQ*', config))
    assert cleavage_gain
    assert deletion_junction
    assert all('*' not in peptide for peptide, _start, _end in stop_gain)


def test_isobaric_normalization_prevents_il_only_novel_peptides():
    exact = PeptideConfig(
        enzyme='No_cut',
        min_length=1,
        max_length=30,
        initiator_methionine='retain',
        isobaric=False,
    )
    isobaric = PeptideConfig(
        enzyme='No_cut',
        min_length=1,
        max_length=30,
        initiator_methionine='retain',
        isobaric=True,
    )
    assert list(iter_alt_specific_peptides('PEPTIDEI', 'PEPTIDEL', exact))
    assert not list(iter_alt_specific_peptides('PEPTIDEI', 'PEPTIDEL', isobaric))


@pytest.mark.parametrize('alias_kind', ['same', 'hardlink', 'symlink'])
def test_build_rejects_annotation_database_aliases(tmp_path, alias_kind):
    source = tmp_path / 'annotation.sqlite'
    make_annotation_sqlite(source)
    original = source.read_bytes()
    if alias_kind == 'same':
        output = source
    else:
        output = tmp_path / f'{alias_kind}.sqlite'
        if alias_kind == 'hardlink':
            os.link(source, output)
        else:
            output.symlink_to(source)

    with pytest.raises(ValueError, match='must be different files'):
        build_peptide_sqlite(source, output, PeptideConfig(), force=True)

    assert source.read_bytes() == original
    connection = sqlite3.connect(source)
    tables = {row[0] for row in connection.execute(
        "SELECT name FROM sqlite_master WHERE type='table'"
    )}
    connection.close()
    assert tables == {'protein_description'}


def test_ensure_rebuild_rejects_annotation_database_as_index(tmp_path):
    source = tmp_path / 'annotation.sqlite'
    make_annotation_sqlite(source)
    with pytest.raises(ValueError, match='must be different files'):
        ensure_peptide_sqlite(source, PeptideConfig(), source, rebuild=True)
    connection = sqlite3.connect(source)
    assert connection.execute('SELECT COUNT(*) FROM protein_description').fetchone()[0] == 3
    connection.close()


def test_builder_cli_rejects_annotation_database_as_output(tmp_path):
    source = tmp_path / 'annotation.sqlite'
    make_annotation_sqlite(source)
    original = source.read_bytes()
    result = subprocess.run(
        [
            sys.executable,
            '-m',
            'precisionprodb.peptideSqlite',
            '-S',
            str(source),
            '-o',
            str(source),
            '--force',
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert 'must be different files' in result.stderr
    assert source.read_bytes() == original


def test_validation_rejects_known_peptide_count_mismatch(index):
    source, output, config, _metadata = index
    connection = sqlite3.connect(output)
    connection.execute('DELETE FROM known_peptide')
    connection.commit()
    connection.close()
    with pytest.raises(ValueError, match='count does not match metadata'):
        validate_peptide_sqlite(output, config, source)


def test_validation_rejects_wrong_user_version(index):
    source, output, config, _metadata = index
    connection = sqlite3.connect(output)
    connection.execute('PRAGMA user_version=2')
    connection.close()
    with pytest.raises(ValueError, match='user_version'):
        validate_peptide_sqlite(output, config, source)


def test_validation_rejects_wrong_table_schema(tmp_path):
    output = tmp_path / 'wrong-schema.sqlite'
    connection = sqlite3.connect(output)
    connection.executescript(
        '''
        PRAGMA user_version=1;
        CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL) WITHOUT ROWID;
        CREATE TABLE known_peptide (wrong_key TEXT PRIMARY KEY) WITHOUT ROWID;
        '''
    )
    connection.close()
    with pytest.raises(ValueError, match='unsupported schema'):
        validate_peptide_sqlite(output, PeptideConfig())


def test_validation_rejects_inconsistent_configuration_hash(index):
    source, output, config, _metadata = index
    update_metadata(output, 'digestion_config_hash', '0' * 64)
    with pytest.raises(ValueError, match='metadata is inconsistent'):
        validate_peptide_sqlite(output, config, source)


def test_validation_wraps_corrupt_sqlite_errors(tmp_path):
    output = tmp_path / 'corrupt.sqlite'
    output.write_bytes(b'not a sqlite database')
    with pytest.raises(ValueError, match='cannot validate peptide SQLite'):
        validate_peptide_sqlite(output, PeptideConfig())


def test_uri_safe_paths_support_sqlite_special_characters(tmp_path):
    source = tmp_path / 'annotation ?#.sqlite'
    output = tmp_path / 'known ?#.peptides.sqlite'
    make_annotation_sqlite(source)
    config = PeptideConfig(min_length=3)
    build_peptide_sqlite(source, output, config)
    validate_peptide_sqlite(output, config, source)
    with KnownPeptideIndex(output) as known:
        assert 'MPEPTIDEK' in known.lookup_many(['MPEPTIDEK'])


def test_enzyme_names_are_canonicalized():
    canonical = PeptideConfig(enzyme='Trypsin')
    lowercase = PeptideConfig(enzyme='trypsin')
    assert lowercase.enzyme == 'Trypsin'
    assert lowercase == canonical
    assert lowercase.config_hash == canonical.config_hash
    assert default_peptide_sqlite_path('annotation.sqlite', lowercase) == (
        default_peptide_sqlite_path('annotation.sqlite', canonical)
    )
    with pytest.raises(ValueError, match='unknown enzyme'):
        PeptideConfig(enzyme='not-an-enzyme')


def test_legacy_lowercase_enzyme_metadata_remains_compatible(index):
    source, output, config, _metadata = index
    legacy_payload = json.loads(config.as_json())
    legacy_payload['enzyme'] = 'trypsin'
    legacy_json = json.dumps(legacy_payload, sort_keys=True, separators=(',', ':'))
    legacy_hash = hashlib.sha256(legacy_json.encode()).hexdigest()
    update_metadata(output, 'digestion_config_json', legacy_json)
    update_metadata(output, 'digestion_config_hash', legacy_hash)

    metadata = validate_peptide_sqlite(output, config, source)
    assert metadata['digestion_config_hash'] == legacy_hash
