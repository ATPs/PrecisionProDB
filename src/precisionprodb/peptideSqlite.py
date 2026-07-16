"""Canonical peptide SQLite indexes used by PrecisionProDB peptide mode."""

import argparse
from dataclasses import asdict, dataclass
import datetime
import hashlib
import json
from multiprocessing import Pool
import os
from pathlib import Path
import sqlite3
import sys

if __package__:
    from .digestion import iter_digest_sequence_occurrences, read_fasta_records
else:
    from digestion import iter_digest_sequence_occurrences, read_fasta_records


SCHEMA_VERSION = 1
DIGESTION_SEMANTICS_VERSION = 1


@dataclass(frozen=True)
class PeptideConfig:
    """The complete digestion profile represented by one peptide index."""

    mode: str = 'enzyme'
    enzyme: str = 'Trypsin'
    missed_cleavages: int = 2
    min_length: int = 7
    max_length: int = 35
    specificity: str = 'full'
    initiator_methionine: str = 'both'
    isobaric: bool = False

    def __post_init__(self):
        if self.mode != 'enzyme':
            raise ValueError("only peptide mode 'enzyme' is supported")
        if self.missed_cleavages < 0:
            raise ValueError('missed_cleavages must be at least 0')
        if self.min_length < 1:
            raise ValueError('min_length must be at least 1')
        if self.max_length < self.min_length:
            raise ValueError('max_length must be at least min_length')
        if self.specificity not in ['full', 'n-specific', 'c-specific', 'semi']:
            raise ValueError("specificity must be 'full', 'n-specific', 'c-specific', or 'semi'")
        if self.initiator_methionine not in ['both', 'remove', 'retain']:
            raise ValueError("initiator_methionine must be 'both', 'remove', or 'retain'")

    @classmethod
    def from_namespace(cls, namespace):
        return cls(
            enzyme=namespace.peptide_enzyme,
            missed_cleavages=namespace.peptide_missed_cleavages,
            min_length=namespace.peptide_min_length,
            max_length=namespace.peptide_max_length,
            specificity=namespace.peptide_specificity,
            initiator_methionine=namespace.peptide_initiator_methionine,
            isobaric=namespace.peptide_isobaric,
        )

    def as_json(self):
        return json.dumps(asdict(self), sort_keys=True, separators=(',', ':'))

    @property
    def config_hash(self):
        return hashlib.sha256(self.as_json().encode()).hexdigest()

    def digest_kwargs(self):
        return {
            'enzyme': self.enzyme,
            'miss_cleavage': self.missed_cleavages,
            'min_len': self.min_length,
            'max_len': self.max_length,
            'specificity': self.specificity,
            'initiator_methionine': self.initiator_methionine,
            'isobaric': self.isobaric,
        }


def default_peptide_sqlite_path(source_sqlite, config):
    """Return a collision-resistant default index path beside the annotation DB."""
    source = Path(source_sqlite)
    stem = source.name[:-7] if source.name.endswith('.sqlite') else source.name
    enzyme = ''.join(char.lower() if char.isalnum() else '-' for char in config.enzyme)
    exactness = 'il' if config.isobaric else 'exact'
    profile = (
        f'{enzyme}_mc{config.missed_cleavages}_len{config.min_length}-{config.max_length}_'
        f'{config.specificity}_im-{config.initiator_methionine}_{exactness}'
    )
    return str(source.with_name(f'{stem}.{profile}.peptides.sqlite'))


def _package_version():
    version_path = Path(__file__).with_name('version')
    try:
        return version_path.read_text().strip()
    except OSError:
        return 'unknown'


def _source_proteome_fingerprint(source_sqlite):
    """Hash a stable protein-id/sequence representation, not SQLite file bytes."""
    source = str(source_sqlite)
    digest = hashlib.sha256()
    protein_count = 0
    uri = f'file:{Path(source).absolute()}?mode=ro'
    connection = sqlite3.connect(uri, uri=True)
    try:
        for protein_id, sequence in connection.execute(
            'SELECT protein_id, AA_seq FROM protein_description ORDER BY protein_id'
        ):
            digest.update(str(protein_id).encode())
            digest.update(b'\0')
            digest.update(str(sequence).encode())
            digest.update(b'\0')
            protein_count += 1
    finally:
        connection.close()
    return digest.hexdigest(), protein_count


def _create_schema(connection):
    connection.executescript(
        '''
        PRAGMA page_size=32768;
        PRAGMA journal_mode=OFF;
        PRAGMA synchronous=OFF;
        PRAGMA temp_store=MEMORY;
        PRAGMA locking_mode=EXCLUSIVE;
        PRAGMA cache_size=-262144;
        CREATE TABLE metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        ) WITHOUT ROWID;
        CREATE TABLE known_peptide (
            peptide_key TEXT PRIMARY KEY
        ) WITHOUT ROWID;
        PRAGMA user_version=1;
        '''
    )


def _iter_peptide_keys(sequence, config):
    for peptide_key, _start, _end in iter_digest_sequence_occurrences(
        sequence, **config.digest_kwargs()
    ):
        yield peptide_key


def _digest_sequence_key_set(payload):
    """Pickle-friendly worker helper for canonical-index construction."""
    sequence, config = payload
    return set(_iter_peptide_keys(sequence, config))


def build_peptide_sqlite(source_sqlite, output_sqlite, peptide_config, threads=1,
                         additional_fastas=None, force=False):
    """Build a complete canonical-peptide index and atomically publish it."""
    if not isinstance(peptide_config, PeptideConfig):
        raise TypeError('peptide_config must be a PeptideConfig')
    source_sqlite = str(source_sqlite)
    output_sqlite = str(output_sqlite)
    if not os.path.exists(source_sqlite):
        raise FileNotFoundError(source_sqlite)
    output_parent = os.path.dirname(os.path.abspath(output_sqlite))
    os.makedirs(output_parent, exist_ok=True)
    if os.path.exists(output_sqlite) and not force:
        raise FileExistsError(f'peptide SQLite already exists: {output_sqlite}')

    temporary = f'{output_sqlite}.building.{os.getpid()}'
    if os.path.exists(temporary):
        os.unlink(temporary)
    source_stat = os.stat(source_sqlite)
    additional_fastas = list(additional_fastas or [])
    additional_metadata = []
    source_hash = hashlib.sha256()
    protein_count = 0
    connection = None
    source_connection = None
    pool = None
    try:
        connection = sqlite3.connect(temporary)
        _create_schema(connection)
        connection.execute('BEGIN IMMEDIATE')
        source_uri = f'file:{Path(source_sqlite).absolute()}?mode=ro'
        # Pool.imap consumes its input generator in a feeder thread.
        source_connection = sqlite3.connect(source_uri, uri=True, check_same_thread=False)
        peptide_batch = set()

        def flush_batch():
            if peptide_batch:
                connection.executemany(
                    'INSERT OR IGNORE INTO known_peptide(peptide_key) VALUES (?)',
                    ((key,) for key in peptide_batch),
                )
                peptide_batch.clear()

        def source_sequences():
            nonlocal protein_count
            for protein_id, sequence in source_connection.execute(
                'SELECT protein_id, AA_seq FROM protein_description ORDER BY protein_id'
            ):
                sequence = str(sequence)
                source_hash.update(str(protein_id).encode())
                source_hash.update(b'\0')
                source_hash.update(sequence.encode())
                source_hash.update(b'\0')
                protein_count += 1
                yield sequence, peptide_config

        worker_count = max(1, min(int(threads), 8))
        if worker_count == 1:
            peptide_sets = map(_digest_sequence_key_set, source_sequences())
        else:
            pool = Pool(worker_count)
            peptide_sets = pool.imap(_digest_sequence_key_set, source_sequences(), chunksize=100)
        for peptide_keys in peptide_sets:
            peptide_batch.update(peptide_keys)
            if len(peptide_batch) >= 50000:
                flush_batch()
        flush_batch()
        if pool is not None:
            pool.close()
            pool.join()
            pool = None

        for fasta_path in additional_fastas:
            fasta_path = str(fasta_path)
            fasta_hash = hashlib.sha256()
            records = 0
            for header, sequence in read_fasta_records(fasta_path):
                fasta_hash.update(header.encode())
                fasta_hash.update(b'\0')
                fasta_hash.update(sequence.encode())
                fasta_hash.update(b'\0')
                records += 1
                peptide_batch.update(_iter_peptide_keys(sequence, peptide_config))
                if len(peptide_batch) >= 50000:
                    flush_batch()
            flush_batch()
            stat = os.stat(fasta_path)
            additional_metadata.append({
                'path': os.path.abspath(fasta_path),
                'size': stat.st_size,
                'mtime_ns': stat.st_mtime_ns,
                'sha256': fasta_hash.hexdigest(),
                'record_count': records,
            })

        known_count = connection.execute('SELECT COUNT(*) FROM known_peptide').fetchone()[0]
        metadata = {
            'schema_version': str(SCHEMA_VERSION),
            'build_complete': '0',
            'digestion_config_json': peptide_config.as_json(),
            'digestion_config_hash': peptide_config.config_hash,
            'digestion_semantics_version': str(DIGESTION_SEMANTICS_VERSION),
            'source_proteome_sha256': source_hash.hexdigest(),
            'source_sqlite_path': os.path.abspath(source_sqlite),
            'source_sqlite_size': str(source_stat.st_size),
            'source_sqlite_mtime_ns': str(source_stat.st_mtime_ns),
            'source_protein_count': str(protein_count),
            'known_peptide_count': str(known_count),
            'additional_fastas_json': json.dumps(additional_metadata, sort_keys=True),
            'precisionprodb_version': _package_version(),
            'created_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        }
        connection.executemany(
            'INSERT INTO metadata(key, value) VALUES (?, ?)', metadata.items()
        )
        connection.execute("UPDATE metadata SET value='1' WHERE key='build_complete'")
        connection.commit()
        connection.close()
        connection = None
        source_connection.close()
        source_connection = None
        os.replace(temporary, output_sqlite)
    except Exception:
        if pool is not None:
            pool.terminate()
            pool.join()
        if connection is not None:
            connection.close()
        if source_connection is not None:
            source_connection.close()
        if os.path.exists(temporary):
            os.unlink(temporary)
        raise
    return validate_peptide_sqlite(output_sqlite, peptide_config, source_sqlite=source_sqlite)


def validate_peptide_sqlite(peptide_sqlite, peptide_config, source_sqlite=None):
    """Validate index structure, profile, completion marker, and source proteome."""
    if not isinstance(peptide_config, PeptideConfig):
        raise TypeError('peptide_config must be a PeptideConfig')
    peptide_sqlite = str(peptide_sqlite)
    if not os.path.exists(peptide_sqlite):
        raise FileNotFoundError(peptide_sqlite)
    uri = f'file:{Path(peptide_sqlite).absolute()}?mode=ro'
    connection = sqlite3.connect(uri, uri=True)
    try:
        tables = {row[0] for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        )}
        if {'metadata', 'known_peptide'} - tables:
            raise ValueError(f'not a peptide SQLite index: {peptide_sqlite}')
        metadata = dict(connection.execute('SELECT key, value FROM metadata'))
    finally:
        connection.close()
    if metadata.get('schema_version') != str(SCHEMA_VERSION):
        raise ValueError('unsupported peptide SQLite schema version')
    if metadata.get('build_complete') != '1':
        raise ValueError('peptide SQLite build is incomplete')
    if metadata.get('digestion_semantics_version') != str(DIGESTION_SEMANTICS_VERSION):
        raise ValueError('peptide SQLite digestion semantics do not match this version')
    if metadata.get('digestion_config_hash') != peptide_config.config_hash:
        raise ValueError('peptide SQLite configuration does not match requested peptide mode')
    if source_sqlite is not None:
        source_hash, source_count = _source_proteome_fingerprint(source_sqlite)
        if source_hash != metadata.get('source_proteome_sha256'):
            raise ValueError('peptide SQLite source proteome does not match annotation SQLite')
        if str(source_count) != metadata.get('source_protein_count'):
            raise ValueError('peptide SQLite source protein count does not match annotation SQLite')
    return metadata


def ensure_peptide_sqlite(source_sqlite, peptide_config, peptide_sqlite=None,
                          rebuild=False, threads=1):
    """Return a valid index, building or atomically rebuilding when requested."""
    peptide_sqlite = str(
        peptide_sqlite or default_peptide_sqlite_path(source_sqlite, peptide_config)
    )
    if rebuild or not os.path.exists(peptide_sqlite):
        build_peptide_sqlite(
            source_sqlite,
            peptide_sqlite,
            peptide_config,
            threads=threads,
            force=os.path.exists(peptide_sqlite),
        )
    else:
        validate_peptide_sqlite(peptide_sqlite, peptide_config, source_sqlite=source_sqlite)
    return peptide_sqlite


class KnownPeptideIndex:
    """Read-only batched membership lookup with a run-local status cache."""

    def __init__(self, peptide_sqlite):
        self.peptide_sqlite = str(peptide_sqlite)
        uri = f'file:{Path(self.peptide_sqlite).absolute()}?mode=ro&immutable=1'
        self.connection = sqlite3.connect(uri, uri=True)
        self.connection.execute('PRAGMA temp_store=MEMORY')
        self.connection.execute('PRAGMA mmap_size=1073741824')
        self.connection.execute(
            'CREATE TEMP TABLE candidate_query (peptide_key TEXT PRIMARY KEY) WITHOUT ROWID'
        )
        self.known_status = {}

    def lookup_many(self, peptide_keys):
        """Return the subset of *peptide_keys* present in the canonical index."""
        unique_keys = {key for key in peptide_keys if key}
        unresolved = [key for key in unique_keys if key not in self.known_status]
        if unresolved:
            self.connection.execute('DELETE FROM candidate_query')
            self.connection.executemany(
                'INSERT OR IGNORE INTO candidate_query(peptide_key) VALUES (?)',
                ((key,) for key in unresolved),
            )
            known = {
                row[0] for row in self.connection.execute(
                    'SELECT c.peptide_key FROM candidate_query AS c '
                    'INNER JOIN known_peptide AS k USING (peptide_key)'
                )
            }
            self.known_status.update({key: key in known for key in unresolved})
        return {key for key in unique_keys if self.known_status[key]}

    def close(self):
        if self.connection is not None:
            self.connection.close()
            self.connection = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


def build_parser():
    parser = argparse.ArgumentParser(
        description='Build an immutable canonical peptide SQLite index from PrecisionProDB annotation SQLite.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('-S', '--sqlite', required=True, help='source annotation SQLite')
    parser.add_argument('-o', '--output', required=True, help='output peptide SQLite')
    parser.add_argument('-t', '--threads', type=int, default=1, help='digestion worker count (maximum 8)')
    parser.add_argument('--enzyme', default='Trypsin')
    parser.add_argument('--missed-cleavages', type=int, default=2)
    parser.add_argument('--min-peptide-length', type=int, default=7)
    parser.add_argument('--max-peptide-length', type=int, default=35)
    parser.add_argument('--specificity', default='full', choices=['full', 'n-specific', 'c-specific', 'semi'])
    parser.add_argument('--initiator-methionine', default='both', choices=['both', 'remove', 'retain'])
    parser.add_argument('--isobaric', action='store_true')
    parser.add_argument('--known-fasta', action='append', default=[], help='additional known FASTA; repeatable')
    parser.add_argument('--force', action='store_true', help='atomically replace an existing output index')
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        config = PeptideConfig(
            enzyme=args.enzyme,
            missed_cleavages=args.missed_cleavages,
            min_length=args.min_peptide_length,
            max_length=args.max_peptide_length,
            specificity=args.specificity,
            initiator_methionine=args.initiator_methionine,
            isobaric=args.isobaric,
        )
        metadata = build_peptide_sqlite(
            args.sqlite,
            args.output,
            config,
            threads=args.threads,
            additional_fastas=args.known_fasta,
            force=args.force,
        )
    except (FileExistsError, FileNotFoundError, ValueError) as exc:
        parser.error(str(exc))
    print(f"built {args.output}: {metadata['known_peptide_count']} known peptides")


if __name__ == '__main__':
    main()
