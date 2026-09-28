import csv
from collections import defaultdict

from Bio import SeqIO

if __package__:
    from .args import add_argument_set
else:
    from args import add_argument_set


def openFile(filename):
    """Open plain or gzipped text files."""
    if filename.endswith('.gz'):
        import gzip
        return gzip.open(filename, 'rt')
    return open(filename)


def _read_fasta_files(files):
    for filename in files.split(','):
        with openFile(filename) as handle:
            yield from SeqIO.parse(handle, 'fasta')


def _normalize_sequence(sequence):
    return str(sequence).strip('*X').upper()


def _write_fasta_entry(handle, uniprot_id, description, sequence, index, count, status='changed'):
    suffix = '' if count == 1 else f'__{index}'
    description = description[len(uniprot_id):].strip()
    header = f'{uniprot_id}{suffix}'
    if description:
        header += f' {description}'
    handle.write(f'>{header}\t{status}\n{sequence}\n')


def extractMutatedUniprot(files_uniprot, files_ref, files_alt, outprefix, length_min=20):
    """Match UniProt proteins to references and report altered sequences.

    Writes ``.uniprot_changed.tsv``, ``.uniprot_changed.fa`` and
    ``.uniprot_all.fa``. Multiple reference IDs with the same sequence are
    retained, and alternatives are sorted before output identifiers are made.
    """
    uniprot_records = []
    for record in _read_fasta_files(files_uniprot):
        uniprot_records.append((record.id, record.description, _normalize_sequence(record.seq)))

    # A protein sequence can be shared by multiple reference identifiers.
    # Preserve every relationship instead of letting FASTA order select one.
    sequence_to_refs = defaultdict(set)
    for record in _read_fasta_files(files_ref):
        sequence = _normalize_sequence(record.seq)
        if len(sequence) >= length_min:
            sequence_to_refs[sequence].add(record.id)

    all_sequences_by_ref = defaultdict(set)
    changed_sequences_by_ref = defaultdict(set)
    for record in _read_fasta_files(files_alt):
        ref_id = record.id.split('__', 1)[0]
        sequence = _normalize_sequence(record.seq)
        all_sequences_by_ref[ref_id].add(sequence)
        if record.description.endswith('\tchanged'):
            changed_sequences_by_ref[ref_id].add(sequence)

    changed_relationships = []
    changed_records = []
    all_records = []
    for uniprot_id, description, uniprot_sequence in uniprot_records:
        ref_ids = sorted(sequence_to_refs.get(uniprot_sequence, ()))
        changed_ref_ids = [ref_id for ref_id in ref_ids if changed_sequences_by_ref.get(ref_id)]
        changed_sequences = sorted({
            sequence
            for ref_id in changed_ref_ids
            for sequence in changed_sequences_by_ref[ref_id]
        })
        if not changed_sequences:
            all_records.append((uniprot_id, description, [uniprot_sequence], False, False, uniprot_sequence))
            continue

        changed_relationships.extend((uniprot_id, ref_id) for ref_id in changed_ref_ids)
        changed_records.append((uniprot_id, description, changed_sequences))

        # The combined file should retain the reference sequence if at least
        # one matching reference ID still has it in the per-genome FASTA.
        reference_is_present = any(
            uniprot_sequence in all_sequences_by_ref.get(ref_id, ())
            for ref_id in ref_ids
        )
        all_records.append((
            uniprot_id, description, changed_sequences, True,
            reference_is_present, uniprot_sequence,
        ))

    with open(outprefix + '.uniprot_changed.tsv', 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['uniprot_id', 'ref_id'])
        writer.writerows(changed_relationships)

    with open(outprefix + '.uniprot_changed.fa', 'w') as changed_handle, open(
        outprefix + '.uniprot_all.fa', 'w'
    ) as all_handle:
        for uniprot_id, description, sequences in changed_records:
            for index, sequence in enumerate(sequences, start=1):
                _write_fasta_entry(
                    changed_handle, uniprot_id, description, sequence, index, len(sequences)
                )

        for (
            uniprot_id, description, sequences, is_changed,
            reference_is_present, uniprot_sequence,
        ) in all_records:
            if not is_changed:
                # Preserve the legacy unchanged header and full UniProt ID.
                all_handle.write(f'>{description}\tunchanged\n{sequences[0]}\n')
                continue
            if reference_is_present:
                uniprot_description = description[len(uniprot_id):].strip()
                header = f'{uniprot_id}__reference'
                if uniprot_description:
                    header += f' {uniprot_description}'
                all_handle.write(f'>{header}\tunchanged\n{uniprot_sequence}\n')
            for index, sequence in enumerate(sequences, start=1):
                _write_fasta_entry(
                    all_handle, uniprot_id, description, sequence, index, len(sequences)
                )


description = '''
Match UniProt proteins in files_uniprot with files_ref.
Output mutated proteins in files_alt.
write three files, outprefix + '.uniprot_changed.tsv'/'.uniprot_changed.fa'/'.uniprot_all.fa'
'''


def build_parser():
    import argparse
    parser = argparse.ArgumentParser(description=description)
    add_argument_set(parser, 'uniprot_matching')
    return parser


def main(argv=None):
    f = build_parser().parse_args(argv)
    extractMutatedUniprot(
        f.files_uniprot,
        f.files_ref,
        f.files_alt,
        f.outprefix,
        length_min=f.length_min,
    )


if __name__ == '__main__':
    main()
