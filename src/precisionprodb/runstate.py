"""Provenance and safe reuse for a PrecisionProDB output prefix."""

import glob
import hashlib
import json
import os
import shutil
import tempfile
import time


FORMAT_VERSION = 2
OUTPUT_SUFFIXES = (
    '.pergeno.protein_all.fa', '.pergeno.protein_changed.fa',
    '.pergeno.aa_mutations.csv', '.pergeno.protein_PEFF.fa',
    '.pergeno.peptide_novel.tsv', '.pergeno.peptide_novel.fa',
    '.pergeno.mutated_protein.fa', '.uniprot_changed.tsv',
    '.uniprot_changed.fa', '.uniprot_all.fa', '.uniprot_PEFF.fa',
    '.file_proteins_input_from_sqlite.fasta',
    '.vcf2mutation.tsv', '.vcf2mutation.tsv.done',
    '.vcf2mutation.tsv.cache.json',
    '.aa_mutations.csv', '.mutated_protein.fa',
    '.vcf2mutation_1.tsv', '.vcf2mutation_2.tsv',
    '_1.pergeno.protein_all.fa', '_1.pergeno.protein_changed.fa',
    '_1.pergeno.aa_mutations.csv',
    '_2.pergeno.protein_all.fa', '_2.pergeno.protein_changed.fa',
    '_2.pergeno.aa_mutations.csv',
)


def file_fingerprint(path):
    """Fingerprint file bytes once; size and mtime are diagnostic only."""
    path = os.path.abspath(path)
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    stat = os.stat(path)
    return {'path': path, 'sha256': digest.hexdigest(), 'size': stat.st_size,
            'mtime_ns': stat.st_mtime_ns}


def mutation_files(value, is_manifest=False):
    if not value:
        return []
    if is_manifest:
        import csv
        with open(value, newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            return [value] + [row['filepath'] for row in reader]
    if os.path.isfile(value):
        return [value]
    if ',' in value and all(os.path.isfile(piece) for piece in value.split(',')):
        return value.split(',')
    matches = sorted(glob.glob(value)) if any(char in value for char in '*?[') else []
    return matches


def _stable_inputs(file_mutations, file_sqlite, file_genome, file_gtf,
                   file_protein, files_uniprot, is_manifest):
    paths = mutation_files(file_mutations, is_manifest=is_manifest)
    if file_sqlite and os.path.isfile(file_sqlite):
        paths.append(file_sqlite)
    else:
        paths.extend(p for p in (file_genome, file_gtf, file_protein) if p)
    paths.extend(p for p in files_uniprot.split(',') if p)
    seen = set()
    result = []
    for path in paths:
        absolute = os.path.abspath(path)
        if absolute not in seen:
            fingerprint = file_fingerprint(path)
            result.append({key: fingerprint[key] for key in ('path', 'sha256', 'size')})
            seen.add(absolute)
    return result


def _generated_paths(outprefix, owned_sqlite_path=None):
    paths = [outprefix + suffix for suffix in OUTPUT_SUFFIXES]
    paths.extend(outprefix + suffix for suffix in ('_temp', '_1_temp', '_2_temp'))
    if owned_sqlite_path:
        paths.append(owned_sqlite_path)
    return [path for path in paths if os.path.lexists(path)]


def sqlite_owned_by_run(outprefix, file_sqlite, file_genome, file_gtf,
                        file_protein, default_mode=False):
    """Distinguish a database this prefix built from an external annotation DB."""
    if not file_sqlite:
        return False
    receipt = outprefix + '.run.json'
    if os.path.isfile(receipt):
        try:
            with open(receipt) as handle:
                previous = json.load(handle)
            settings = previous.get('settings', {})
            if settings.get('sqlite') == os.path.abspath(file_sqlite):
                return bool(settings.get('owned_sqlite'))
        except (OSError, ValueError):
            pass  # RunState.prepare reports the invalid receipt.
    return bool(file_genome and file_gtf and file_protein and
                (not os.path.exists(file_sqlite) or default_mode))


def _write_json(path, value):
    temporary = path + f'.tmp.{os.getpid()}'
    with open(temporary, 'w') as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _temporary_artifacts(outprefix):
    result = {}
    for suffix in ('_temp', '_1_temp', '_2_temp'):
        folder = outprefix + suffix
        if not os.path.isdir(folder):
            continue
        for root, _directories, files in os.walk(folder):
            for name in files:
                path = os.path.join(root, name)
                stat = os.stat(path)
                result[os.path.abspath(path)] = [stat.st_size, stat.st_mtime_ns]
    return result


class RunState:
    def __init__(self, outprefix, inputs, settings, force=False, keep_all=False,
                 owned_sqlite=False):
        self.outprefix = outprefix
        self.path = outprefix + '.run.json'
        self.inputs = inputs
        self.settings = settings
        self.force = force
        self.keep_all = keep_all
        self.owned_sqlite = owned_sqlite
        self.owned_sqlite_path = (settings.get('sqlite') or outprefix + '.sqlite') if owned_sqlite else None
        self.current = {
            'format_version': FORMAT_VERSION,
            'inputs': inputs,
            'settings': settings,
        }

    @classmethod
    def from_inputs(cls, outprefix, file_mutations, file_sqlite, file_genome,
                    file_gtf, file_protein, files_uniprot, is_manifest,
                    settings, force=False, keep_all=False, owned_sqlite=False):
        if owned_sqlite:
            file_sqlite = ''
        inputs = _stable_inputs(file_mutations, file_sqlite, file_genome,
                                file_gtf, file_protein, files_uniprot, is_manifest)
        return cls(outprefix, inputs, settings, force=force, keep_all=keep_all,
                   owned_sqlite=owned_sqlite)

    def prepare(self):
        previous = None
        if os.path.exists(self.path):
            try:
                with open(self.path) as handle:
                    previous = json.load(handle)
            except (OSError, ValueError) as exc:
                raise ValueError(f'unreadable run receipt {self.path}: {exc}') from exc
        generated = _generated_paths(self.outprefix, self.owned_sqlite_path)
        if self.force:
            if generated or previous:
                if self.keep_all:
                    archive_base = self.outprefix + '.archive'
                    os.makedirs(archive_base, exist_ok=True)
                    archive = tempfile.mkdtemp(
                        prefix=time.strftime('%Y%m%d-%H%M%S-'), dir=archive_base
                    )
                for path in generated + ([self.path] if previous else []):
                    if self.keep_all:
                        shutil.move(path, os.path.join(archive, os.path.basename(path)))
                    elif os.path.isdir(path) and not os.path.islink(path):
                        shutil.rmtree(path)
                    else:
                        os.unlink(path)
        elif previous is None and generated:
            raise ValueError(f'unverified generated files exist for {self.outprefix}; use a fresh prefix or --force')
        elif previous is not None:
            for key, value in self.current.items():
                if previous.get(key) != value:
                    raise ValueError(f'run inputs/settings changed for {self.outprefix}; use a fresh prefix or --force')
            if previous.get('status') != 'complete':
                raise ValueError(f'previous run for {self.outprefix} was incomplete; use --force to rebuild')
            for path, expected in previous.get('outputs', {}).items():
                if not os.path.exists(path) or file_fingerprint(path)['sha256'] != expected:
                    raise ValueError(f'generated output changed or is missing: {path}; use --force')
            if _temporary_artifacts(self.outprefix) != previous.get('temporary_artifacts', {}):
                raise ValueError(f'intermediate files changed for {self.outprefix}; use --force')
            return True
        self.record('running')
        return False

    def record(self, status, **details):
        data = dict(self.current, status=status, updated_at=time.time(), **details)
        _write_json(self.path, data)

    def complete(self):
        required = [
            self.outprefix + '.pergeno.aa_mutations.csv',
            self.outprefix + '.pergeno.protein_all.fa',
            self.outprefix + '.pergeno.protein_changed.fa',
        ]
        if self.settings.get('peff'):
            required.append(self.outprefix + '.pergeno.protein_PEFF.fa')
        if self.settings.get('peptide'):
            required.extend([
                self.outprefix + '.pergeno.peptide_novel.tsv',
                self.outprefix + '.pergeno.peptide_novel.fa',
            ])
        if (self.settings.get('download') == 'UNIPROT' and
                mutation_files(self.settings.get('mutations', ''))):
            required.extend([
                self.outprefix + '.uniprot_changed.tsv',
                self.outprefix + '.uniprot_changed.fa',
                self.outprefix + '.uniprot_all.fa',
            ])
            if self.settings.get('peff'):
                required.append(self.outprefix + '.uniprot_PEFF.fa')
        missing = [path for path in required if not os.path.isfile(path)]
        if missing:
            raise RuntimeError('run did not produce required outputs: ' + ', '.join(missing))
        outputs = {}
        counts = {}
        for path in _generated_paths(self.outprefix, self.owned_sqlite_path):
            if not os.path.isfile(path):
                continue
            outputs[os.path.abspath(path)] = file_fingerprint(path)['sha256']
            if path.endswith(('.fa', '.fasta')):
                with open(path) as handle:
                    counts[os.path.basename(path)] = sum(line.startswith('>') for line in handle)
            elif path.endswith(('.csv', '.tsv')):
                with open(path) as handle:
                    counts[os.path.basename(path)] = max(0, sum(1 for _ in handle) - 1)
        self.record('complete', outputs=outputs, counts=counts,
                    temporary_artifacts=_temporary_artifacts(self.outprefix))

    def fail(self, exc):
        self.record('failed', error=f'{type(exc).__name__}: {exc}')


class StageCache:
    """Validate an individual converter result before trusting its .done file."""

    def __init__(self, output, input_files, settings, force=False):
        self.output = output
        self.done = output + '.done'
        self.receipt = output + '.cache.json'
        self.force = force
        self.current = {
            'format_version': FORMAT_VERSION,
            'inputs': [
                {key: fp[key] for key in ('path', 'sha256', 'size')}
                for fp in (file_fingerprint(path) for path in input_files)
            ],
            'settings': settings,
        }

    def prepare(self):
        paths = [p for p in (self.output, self.done, self.receipt)
                 if os.path.lexists(p)]
        if self.force and paths:
            archive_base = self.output + '.archive'
            os.makedirs(archive_base, exist_ok=True)
            archive = tempfile.mkdtemp(
                prefix=time.strftime('%Y%m%d-%H%M%S-'), dir=archive_base
            )
            for path in paths:
                shutil.move(path, os.path.join(archive, os.path.basename(path)))
            return False
        if not paths:
            return False
        if len(paths) != 3:
            raise ValueError(f'incomplete cached converter output {self.output}; use --force')
        try:
            with open(self.receipt) as handle:
                previous = json.load(handle)
        except (OSError, ValueError) as exc:
            raise ValueError(f'invalid converter receipt {self.receipt}; use --force') from exc
        if any(previous.get(k) != v for k, v in self.current.items()):
            raise ValueError(f'converter inputs/settings changed for {self.output}; use --force')
        for path in (self.output, self.done):
            if file_fingerprint(path)['sha256'] != previous.get('artifacts', {}).get(os.path.abspath(path)):
                raise ValueError(f'cached converter artifact changed: {path}; use --force')
        return True

    def complete(self):
        artifacts = {
            os.path.abspath(path): file_fingerprint(path)['sha256']
            for path in (self.output, self.done)
        }
        _write_json(self.receipt, dict(self.current, artifacts=artifacts))
