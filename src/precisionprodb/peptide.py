"""Novel altered-protein peptide extraction and deterministic output merging."""

import csv
from multiprocessing import Pool
import os

if __package__:
    from .digestion import iter_digest_sequence_occurrences
else:
    from digestion import iter_digest_sequence_occurrences


PART_COLUMNS = [
    'peptide_sequence', 'protein_id_fasta', 'protein_id', 'individual', 'seqname',
    'peptide_start_alt', 'peptide_end_alt', 'variant_AA', 'insertion_AA',
    'deletion_AA', 'frameChange', 'stopGain', 'AA_stopGain', 'stopLoss',
    'stopLoss_pos', 'enzyme', 'max_missed_cleavages', 'min_peptide_length',
    'max_peptide_length', 'specificity', 'initiator_methionine', 'isobaric',
    'peptide_config_hash',
]
FINAL_COLUMNS = ['peptide_id'] + PART_COLUMNS
_PEPTIDE_WORKER_CONFIG = None


def _digest_occurrences(sequence, config):
    return list(iter_digest_sequence_occurrences(sequence, **config.digest_kwargs()))


def iter_alt_specific_peptides(ref_sequence, alt_sequence, config):
    """Yield altered peptide occurrences absent from the same reference protein."""
    reference_keys = {peptide for peptide, _start, _end in _digest_occurrences(ref_sequence, config)}
    for peptide, start, end in _digest_occurrences(alt_sequence, config):
        if peptide not in reference_keys:
            yield peptide, start, end


def _value(row, name, default=''):
    value = row.get(name, default)
    if value is None:
        return ''
    try:
        if value != value:  # pandas NA/NaN
            return ''
    except (TypeError, ValueError):
        pass
    return value


def _part_record(row, peptide, start, end, config):
    return {
        'peptide_sequence': peptide,
        'protein_id_fasta': _value(row, 'protein_id_fasta'),
        'protein_id': _value(row, 'protein_id'),
        'individual': _value(row, 'individual'),
        'seqname': _value(row, 'seqname'),
        'peptide_start_alt': start,
        'peptide_end_alt': end,
        'variant_AA': _value(row, 'variant_AA'),
        'insertion_AA': _value(row, 'insertion_AA'),
        'deletion_AA': _value(row, 'deletion_AA'),
        'frameChange': _value(row, 'frameChange'),
        'stopGain': _value(row, 'stopGain'),
        'AA_stopGain': _value(row, 'AA_stopGain'),
        'stopLoss': _value(row, 'stopLoss'),
        'stopLoss_pos': _value(row, 'stopLoss_pos'),
        'enzyme': config.enzyme,
        'max_missed_cleavages': config.missed_cleavages,
        'min_peptide_length': config.min_length,
        'max_peptide_length': config.max_length,
        'specificity': config.specificity,
        'initiator_methionine': config.initiator_methionine,
        'isobaric': config.isobaric,
        'peptide_config_hash': config.config_hash,
    }


def _candidate_records_for_group(rows, config):
    """Digest one canonical-protein group without touching the SQLite index."""
    if not rows:
        return []
    reference_sequence = str(rows[0]['AA_seq'])
    reference_keys = {
        peptide for peptide, _start, _end in _digest_occurrences(reference_sequence, config)
    }
    alt_cache = {}
    records = []
    for row in rows:
        alt_sequence = _value(row, 'new_AA')
        if not alt_sequence or alt_sequence == reference_sequence:
            continue
        if alt_sequence not in alt_cache:
            alt_cache[alt_sequence] = [
                (peptide, start, end)
                for peptide, start, end in _digest_occurrences(alt_sequence, config)
                if peptide not in reference_keys
            ]
        for peptide, start, end in alt_cache[alt_sequence]:
            records.append(_part_record(row, peptide, start, end, config))
    return records


def _init_peptide_worker(config):
    global _PEPTIDE_WORKER_CONFIG
    _PEPTIDE_WORKER_CONFIG = config


def _candidate_records_worker(rows):
    return _candidate_records_for_group(rows, _PEPTIDE_WORKER_CONFIG)


def write_novel_peptide_part(df_changed, known_index, output_tsv, config,
                             buffer_size=20000, threads=1):
    """Write novel peptide mappings for one prepared changed-protein dataframe."""
    output_dir = os.path.dirname(os.path.abspath(output_tsv))
    os.makedirs(output_dir, exist_ok=True)
    with open(output_tsv, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=PART_COLUMNS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        if df_changed is None or df_changed.shape[0] == 0:
            return output_tsv
        if 'protein_id' not in df_changed.columns:
            raise ValueError('changed-protein rows must contain protein_id')
        required = {'AA_seq', 'new_AA', 'protein_id_fasta'}
        missing = required - set(df_changed.columns)
        if missing:
            raise ValueError(f'changed-protein rows missing required columns: {sorted(missing)}')

        pending = []

        def flush_pending():
            if not pending:
                return
            known = known_index.lookup_many(item['_peptide_key'] for item in pending)
            for item in pending:
                if item['_peptide_key'] not in known:
                    item.pop('_peptide_key')
                    writer.writerow(item)
            pending.clear()

        groups = (
            group.to_dict(orient='records')
            for _protein_id, group in df_changed.groupby('protein_id', sort=False)
        )
        worker_count = max(1, min(int(threads), 8))
        pool = None
        try:
            if worker_count == 1:
                candidate_groups = (
                    _candidate_records_for_group(rows, config) for rows in groups
                )
            else:
                pool = Pool(worker_count, initializer=_init_peptide_worker, initargs=(config,))
                candidate_groups = pool.imap(_candidate_records_worker, groups, chunksize=32)
            for candidate_records in candidate_groups:
                for record in candidate_records:
                    record['_peptide_key'] = record['peptide_sequence']
                    pending.append(record)
                    if len(pending) >= buffer_size:
                        flush_pending()
            flush_pending()
            if pool is not None:
                pool.close()
                pool.join()
                pool = None
        except Exception:
            if pool is not None:
                pool.terminate()
                pool.join()
            raise
    return output_tsv


def merge_novel_peptide_parts(files, outprefix):
    """Merge TSV parts into deterministic final TSV and nonredundant FASTA."""
    files = [str(path) for path in files if path and os.path.exists(path)]
    final_tsv = outprefix + '.pergeno.peptide_novel.tsv'
    final_fasta = outprefix + '.pergeno.peptide_novel.fa'
    output_dir = os.path.dirname(os.path.abspath(final_tsv))
    os.makedirs(output_dir, exist_ok=True)

    mapping_counts = {}
    for path in files:
        with open(path, newline='') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                peptide = row['peptide_sequence']
                mapping_counts[peptide] = mapping_counts.get(peptide, 0) + 1
    peptide_ids = {
        peptide: f'PPDBpep_{index:06d}'
        for index, peptide in enumerate(sorted(mapping_counts), start=1)
    }

    with open(final_tsv, 'w', newline='') as tsv_handle, open(final_fasta, 'w') as fasta_handle:
        writer = csv.DictWriter(tsv_handle, fieldnames=FINAL_COLUMNS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for peptide in sorted(mapping_counts):
            fasta_handle.write(
                f'>{peptide_ids[peptide]}|n_mappings={mapping_counts[peptide]}\n{peptide}\n'
            )
        for path in files:
            with open(path, newline='') as handle:
                for row in csv.DictReader(handle, delimiter='\t'):
                    row['peptide_id'] = peptide_ids[row['peptide_sequence']]
                    writer.writerow(row)
    return final_tsv, final_fasta
