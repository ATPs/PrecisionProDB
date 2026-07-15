import argparse
import gzip
import sys


ALL_RESIDUES = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
DEFAULT_ENZYME = 'Trypsin'
ENZYME_PRESETS = {
    'Cut_everywhere': {'pos': 'n', 'sites': ALL_RESIDUES, 'no': ''},
    'Trypsin': {'pos': 'c', 'sites': 'KR', 'no': 'P'},
    'Trypsin/P': {'pos': 'c', 'sites': 'KR', 'no': ''},
    'Lys_C': {'pos': 'c', 'sites': 'K', 'no': 'P'},
    'Lys_N': {'pos': 'n', 'sites': 'K', 'no': ''},
    'Arg_C': {'pos': 'c', 'sites': 'R', 'no': 'P'},
    'Asp_N': {'pos': 'n', 'sites': 'DN', 'no': ''},
    'CNBr': {'pos': 'c', 'sites': 'M', 'no': ''},
    'Asp-N_ambic': {'pos': 'c', 'sites': 'DE', 'no': ''},
    'PepsinA': {'pos': 'c', 'sites': 'FL', 'no': ''},
    'Chymotrypsin': {'pos': 'c', 'sites': 'FWYL', 'no': 'P'},
    'No_cut': {'pos': 'c', 'sites': '@', 'no': '@'},
}
ENZYME_LOOKUP = {name.lower(): name for name in ENZYME_PRESETS}
DESCRIPTION = """Digest protein sequences from FASTA input or a direct sequence string.

Input modes:
  1. FASTA file(s):
       python src/precisionprodb/digestion.py proteins.fa proteins2.fa.gz
  2. One sequence:
       python src/precisionprodb/digestion.py --sequence AKRPQK

Cleavage rules:
  - Default behavior is the same as --enzyme Trypsin --missed-cleavages 2.
  - Use --enzyme for Comet-style presets.
  - Use -c/-a/-p to override the preset manually.

I/L handling in this standalone CLI:
  - Default: keep I and L distinct.
  - --isobaric: normalize I -> L before digestion.

Initiator methionine handling:
  - --initiator-methionine controls whether an N-terminal M is retained,
    removed, or both forms are digested (default: both).

Stop-codon handling:
  - Internal "*" is always treated as a hard breakpoint.
  - Output peptides never include "*".
  - Internal "*" still counts in TSV peptide coordinates.

Output behavior:
  - TSV always includes protein_id, peptide, start, and end.
  - protein_id is the first whitespace-delimited token from the FASTA header.
  - peptide and fasta outputs are de-duplicated by peptide sequence.
"""
EPILOG = """Examples:
  Single-sequence digest:
    python src/precisionprodb/digestion.py --sequence AKRPQK -l 2

  FASTA digestion with missed cleavages:
    python src/precisionprodb/digestion.py proteins.fa --enzyme Trypsin --missed-cleavages 2 -l 6 -M 40

  Gzipped FASTA input:
    python src/precisionprodb/digestion.py proteins.fa.gz --output-format peptide

  Remove the initiator methionine before digestion:
    python src/precisionprodb/digestion.py --sequence MAKRPQK --initiator-methionine remove

  Output-format examples:
    python src/precisionprodb/digestion.py proteins.fa --output-format tsv
      sp|P12345\tAK\t1\t2
      sp|P12345\tRPQK\t3\t6

    python src/precisionprodb/digestion.py --sequence 'AK*RPQK' --output-format tsv -l 1
      sequence\tAK\t1\t2
      sequence\tRPQK\t4\t7

    python src/precisionprodb/digestion.py proteins.fa --output-format peptide
      AK
      RPQK

    python src/precisionprodb/digestion.py proteins.fa --output-format fasta
      >sp|P12345_1
      AK
      >sp|P12345_2
      RPQK

  Header-template examples (only used with --output-format fasta):
    --header-template "{protein_id}_{index}"     -> >sp|P12345_1
    --header-template "{protein_id}|pep{index}" -> >sp|P12345|pep1
    --header-template "pep{index}"              -> >pep1

  Enzyme preset examples:
    python src/precisionprodb/digestion.py proteins.fa --enzyme Trypsin
    python src/precisionprodb/digestion.py proteins.fa --enzyme No_cut
    python src/precisionprodb/digestion.py proteins.fa --enzyme Cut_everywhere
    python src/precisionprodb/digestion.py proteins.fa --enzyme Trypsin -a ""

Comet-style enzyme presets for --enzyme NAME:
  0.  Cut_everywhere         0      -           -
  1.  Trypsin                1      KR          P
  2.  Trypsin/P              1      KR          -
  3.  Lys_C                  1      K           P
  4.  Lys_N                  0      K           -
  5.  Arg_C                  1      R           P
  6.  Asp_N                  0      DN          -
  7.  CNBr                   1      M           -
  8.  Asp-N_ambic            1      DE          -
  9.  PepsinA                1      FL          -
  10. Chymotrypsin           1      FWYL        P
  11. No_cut                 1      @           @

Notes:
  - In the preset table, 1 means c-terminal cleavage and 0 means n-terminal cleavage.
  - Cut_everywhere emits all contiguous peptides within the requested length range.
  - Cut_everywhere behaves like unlimited missed cleavages and ignores --missed-cleavages.
  - No_cut emits the full sequence as one peptide.
  - Manual -c/-a/-p values override --enzyme when both are present.
  - Internal "*" always splits the protein before digestion and is never kept in output peptides.
  - TSV coordinates are 1-based inclusive and internal "*" counts in the parent protein position.
  - peptide and fasta outputs keep only the first occurrence of each peptide sequence.
  - --initiator-methionine both emits retained and N-terminal-M-removed forms;
    shared peptide occurrences are emitted once.
  - This PrecisionProDB version is pure Python and does not use fast_digest/Cython.
"""

def open_text_file(filename):
    """Open a plain-text or gzipped text file."""
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    return open(filename)


def read_fasta_records(file_path):
    """
    Read a FASTA file and yield ``(header, sequence)`` pairs.

    The header includes ``>`` and the sequence has leading/trailing ``*`` removed.
    """
    with open_text_file(file_path) as file_handle:
        header = ''
        sequence_parts = []
        for line in file_handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if sequence_parts:
                    yield header, ''.join(sequence_parts).strip('*')
                    sequence_parts = []
                header = line
            else:
                sequence_parts.append(line)
        if sequence_parts:
            yield header, ''.join(sequence_parts).strip('*')


def normalize_sequence(sequence, isobaric=False):
    """Normalize a protein sequence before digestion."""
    sequence = str(sequence).upper().strip('*')
    if isobaric:
        sequence = sequence.replace('I', 'L')
    return sequence


def iter_stop_segments(sequence):
    """Yield ``(start_0based, segment_sequence)`` pairs split on internal stop markers."""
    segment_start = 0
    for index, residue in enumerate(sequence):
        if residue != '*':
            continue
        if segment_start < index:
            yield segment_start, sequence[segment_start:index]
        segment_start = index + 1
    if segment_start < len(sequence):
        yield segment_start, sequence[segment_start:]


def iter_initiator_methionine_variants(sequence, initiator_methionine='both'):
    """Yield ``(start_offset, sequence)`` variants for initiator methionine handling."""
    if initiator_methionine not in ['both', 'remove', 'retain']:
        raise ValueError("initiator_methionine must be 'both', 'remove', or 'retain'")
    if not sequence.startswith('M') or initiator_methionine == 'retain':
        yield 0, sequence
        return
    if initiator_methionine == 'remove':
        yield 1, sequence[1:]
        return
    yield 0, sequence
    yield 1, sequence[1:]


def _is_cut_everywhere_rule(sites, pos, no):
    """Return True when the cleavage rule means every position can be a cut site."""
    return sites == ALL_RESIDUES and pos == 'n' and no == ''


def resolve_enzyme_name(name):
    """Return the canonical preset name for an enzyme."""
    if not name:
        return None
    resolved = ENZYME_LOOKUP.get(str(name).lower())
    if resolved is None:
        raise ValueError(
            'unknown enzyme {!r}. Choose from: {}'.format(
                name,
                ', '.join(ENZYME_PRESETS.keys()),
            )
        )
    return resolved


def resolve_cleavage_rule(enzyme=DEFAULT_ENZYME, cleavage_sites=None, anti_cleavage_sites=None, cleavage_position=None):
    """Resolve preset and manual enzyme settings into one cleavage rule."""
    enzyme_name = resolve_enzyme_name(enzyme) or DEFAULT_ENZYME
    preset = ENZYME_PRESETS[enzyme_name]
    sites = preset['sites'] if cleavage_sites is None else cleavage_sites
    anti_sites = preset['no'] if anti_cleavage_sites is None else anti_cleavage_sites
    position = preset['pos'] if cleavage_position is None else cleavage_position
    if position not in ['c', 'n']:
        raise ValueError("cleavage_position must be 'c' or 'n'")
    return {
        'enzyme': enzyme_name,
        'sites': sites,
        'no': anti_sites,
        'pos': position,
    }


def _unique_positions(positions):
    unique_positions = []
    for position in positions:
        if not unique_positions or unique_positions[-1] != position:
            unique_positions.append(position)
    return unique_positions


def _cleavage_boundaries(protein, sites='KR', pos='c', no='P'):
    if pos not in ['c', 'n']:
        raise ValueError("pos must be 'c' or 'n'")

    boundaries = [0]
    protein_length = len(protein)
    if protein_length == 0:
        return boundaries

    if pos == 'c':
        for index, residue in enumerate(protein):
            if residue not in sites:
                continue
            if index + 1 < protein_length and protein[index + 1] in no:
                continue
            boundaries.append(index + 1)
    else:
        for index, residue in enumerate(protein):
            if residue not in sites:
                continue
            if index + 1 < protein_length and protein[index + 1] in no:
                continue
            boundaries.append(index)

    boundaries.append(protein_length)
    return _unique_positions(boundaries)


def _segment_digest_occurrences(protein, sites='KR', pos='c', no='P', min_len=0, max_len=None):
    """Return peptide occurrences as ``(peptide, start_0based, end_0based_exclusive)``."""
    occurrences = []
    boundaries = _cleavage_boundaries(protein, sites=sites, pos=pos, no=no)
    for start, end in zip(boundaries, boundaries[1:]):
        peptide = protein[start:end]
        peptide_length = len(peptide)
        if not peptide or peptide_length < min_len:
            continue
        if max_len is not None and peptide_length > max_len:
            continue
        occurrences.append((peptide, start, end))
    return occurrences


def _segment_cut_everywhere_occurrences(protein, min_len=0, max_len=None):
    """Return all contiguous peptide occurrences for the requested length range."""
    occurrences = []
    protein_length = len(protein)
    min_len = max(1, min_len)
    if max_len is None:
        max_len = protein_length
    max_len = min(max_len, protein_length)
    for start in range(protein_length):
        max_end = min(protein_length, start + max_len)
        min_end = start + min_len
        if min_end > max_end:
            continue
        for end in range(min_end, max_end + 1):
            occurrences.append((protein[start:end], start, end))
    return occurrences


def _segment_digest_with_missed_cleavages_occurrences(protein, sites='KR', pos='c', no='P', miss_cleavage=2, peplen_min=6, peplen_max=40):
    """Return missed-cleavage occurrences as ``(peptide, start_0based, end_0based_exclusive)``."""
    peptides_cut_all = _segment_digest_occurrences(protein, sites=sites, pos=pos, no=no, min_len=0)
    occurrences = []
    max_missed = max(0, int(miss_cleavage))
    for missed_count in range(max_missed + 1):
        for start_index in range(len(peptides_cut_all) - missed_count):
            segment_occurrences = peptides_cut_all[start_index:start_index + missed_count + 1]
            peptide = ''.join(item[0] for item in segment_occurrences)
            peptide_length = len(peptide)
            if peplen_min <= peptide_length <= peplen_max:
                occurrences.append(
                    (
                        peptide,
                        segment_occurrences[0][1],
                        segment_occurrences[-1][2],
                    )
                )
    return occurrences


def digest(protein, sites='KR', pos='c', no='P', min_len=0, max_len=None):
    """Return non-overlapping peptides after cleavage."""
    peptides = []
    for _segment_start, segment in iter_stop_segments(protein):
        if _is_cut_everywhere_rule(sites, pos, no):
            peptides.extend(
                peptide
                for peptide, _start, _end in _segment_cut_everywhere_occurrences(
                    segment,
                    min_len=min_len,
                    max_len=max_len,
                )
            )
            continue
        peptides.extend(
            peptide
            for peptide, _start, _end in _segment_digest_occurrences(
                segment,
                sites=sites,
                pos=pos,
                no=no,
                min_len=min_len,
                max_len=max_len,
            )
        )
    return peptides


def digest_with_missed_cleavages(protein, sites='KR', pos='c', no='P', miss_cleavage=2, peplen_min=6, peplen_max=40):
    """Return peptides including missed-cleavage combinations."""
    peptides = []
    for _segment_start, segment in iter_stop_segments(protein):
        if _is_cut_everywhere_rule(sites, pos, no):
            peptides.extend(
                peptide
                for peptide, _start, _end in _segment_cut_everywhere_occurrences(
                    segment,
                    min_len=peplen_min,
                    max_len=peplen_max,
                )
            )
            continue
        peptides.extend(
            peptide
            for peptide, _start, _end in _segment_digest_with_missed_cleavages_occurrences(
                segment,
                sites=sites,
                pos=pos,
                no=no,
                miss_cleavage=miss_cleavage,
                peplen_min=peplen_min,
                peplen_max=peplen_max,
            )
        )
    return peptides


def digest_sequence_occurrences(sequence, enzyme=DEFAULT_ENZYME, cleavage_sites=None, anti_cleavage_sites=None, cleavage_position=None, miss_cleavage=2, min_len=6, max_len=40, isobaric=False, initiator_methionine='both'):
    """Digest one protein sequence with missed cleavages and return peptide occurrences."""
    rule = resolve_cleavage_rule(
        enzyme=enzyme,
        cleavage_sites=cleavage_sites,
        anti_cleavage_sites=anti_cleavage_sites,
        cleavage_position=cleavage_position,
    )
    normalized_sequence = normalize_sequence(sequence, isobaric=isobaric)
    peptide_occurrences = []
    seen_occurrences = set()
    for sequence_offset, sequence_variant in iter_initiator_methionine_variants(
        normalized_sequence,
        initiator_methionine=initiator_methionine,
    ):
        for segment_offset, segment in iter_stop_segments(sequence_variant):
            if _is_cut_everywhere_rule(rule['sites'], rule['pos'], rule['no']):
                segment_occurrences = _segment_cut_everywhere_occurrences(
                    segment,
                    min_len=min_len,
                    max_len=max_len,
                )
            else:
                segment_occurrences = _segment_digest_with_missed_cleavages_occurrences(
                    segment,
                    sites=rule['sites'],
                    pos=rule['pos'],
                    no=rule['no'],
                    miss_cleavage=miss_cleavage,
                    peplen_min=min_len,
                    peplen_max=max_len,
                )
            for peptide, start, end in segment_occurrences:
                occurrence = (
                    peptide,
                    sequence_offset + segment_offset + start + 1,
                    sequence_offset + segment_offset + end,
                )
                if occurrence not in seen_occurrences:
                    seen_occurrences.add(occurrence)
                    peptide_occurrences.append(occurrence)
    return peptide_occurrences


def digest_sequence(sequence, enzyme=DEFAULT_ENZYME, cleavage_sites=None, anti_cleavage_sites=None, cleavage_position=None, miss_cleavage=2, min_len=6, max_len=40, isobaric=False, initiator_methionine='both'):
    """Digest one protein sequence with missed cleavages."""
    return [
        peptide
        for peptide, _start, _end in digest_sequence_occurrences(
            sequence,
            enzyme=enzyme,
            cleavage_sites=cleavage_sites,
            anti_cleavage_sites=anti_cleavage_sites,
            cleavage_position=cleavage_position,
            miss_cleavage=miss_cleavage,
            min_len=min_len,
            max_len=max_len,
            isobaric=isobaric,
            initiator_methionine=initiator_methionine,
        )
    ]


def iter_digest_records(file_path, enzyme=DEFAULT_ENZYME, cleavage_sites=None, anti_cleavage_sites=None, cleavage_position=None, miss_cleavage=2, min_len=6, max_len=40, isobaric=False, initiator_methionine='both'):
    """Yield ``(header, peptide_occurrences)`` pairs for a FASTA file."""
    for header, sequence in read_fasta_records(file_path):
        yield header, digest_sequence_occurrences(
            sequence,
            enzyme=enzyme,
            cleavage_sites=cleavage_sites,
            anti_cleavage_sites=anti_cleavage_sites,
            cleavage_position=cleavage_position,
            miss_cleavage=miss_cleavage,
            min_len=min_len,
            max_len=max_len,
            isobaric=isobaric,
            initiator_methionine=initiator_methionine,
        )


def format_tsv_protein_id(header):
    """Format the TSV protein ID from a FASTA header or synthetic sequence name."""
    header_text = header[1:] if header.startswith('>') else header
    if not header_text:
        return header_text
    return header_text.split(None, 1)[0]


def format_fasta_header(header, index, template='{protein_id}_{index}'):
    """Format one peptide FASTA header from a source protein header."""
    header_text = header[1:] if header.startswith('>') else header
    protein_id = format_tsv_protein_id(header)
    try:
        return template.format(protein_id=protein_id, index=index, header=header_text)
    except KeyError as exc:
        raise ValueError('unsupported header template field: {}'.format(exc.args[0]))


def write_digest_results(results, handle, output_format='tsv', header_template='{protein_id}_{index}'):
    """Write digested peptides in TSV, peptide-only, or FASTA format."""
    if output_format not in ['tsv', 'peptide', 'fasta']:
        raise ValueError("output_format must be 'tsv', 'peptide', or 'fasta'")

    seen_peptides = set()
    header_indices = {}
    for header, peptide_occurrences in results:
        if output_format == 'tsv':
            protein_id = format_tsv_protein_id(header)
            for peptide, start, end in peptide_occurrences:
                handle.write('{}\t{}\t{}\t{}\n'.format(protein_id, peptide, start, end))
            continue
        if output_format == 'peptide':
            for peptide, _start, _end in peptide_occurrences:
                if peptide in seen_peptides:
                    continue
                seen_peptides.add(peptide)
                handle.write(peptide + '\n')
            continue
        for peptide, _start, _end in peptide_occurrences:
            if peptide in seen_peptides:
                continue
            seen_peptides.add(peptide)
            index = header_indices.get(header, 0) + 1
            header_indices[header] = index
            fasta_header = format_fasta_header(header, index, template=header_template)
            handle.write('>{}\n{}\n'.format(fasta_header, peptide))


def build_parser():
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        epilog=EPILOG,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        'fasta',
        metavar='*.fasta|*.fa|*.fasta.gz|*.fa.gz|*.txt|*.txt.gz',
        nargs='*',
        help='FASTA file(s) containing protein sequences to digest',
    )
    parser.add_argument(
        '--sequence',
        default='',
        help='Digest a single protein sequence directly instead of reading FASTA input',
    )
    parser.add_argument(
        '--enzyme',
        default='',
        help='Comet-style enzyme preset name.\n'
             'Examples: Trypsin, Trypsin/P, Lys_N, No_cut, Cut_everywhere',
    )
    parser.add_argument(
        '--cleavage_sites',
        '-c',
        dest='cleavage_sites',
        default=None,
        help='Manual cleavage residues. Overrides --enzyme.\n'
             'Example: -c KR',
    )
    parser.add_argument(
        '--anti_cleavage_sites',
        '-a',
        dest='anti_cleavage_sites',
        default=None,
        help='Manual anti-cleavage residues. Overrides --enzyme.\n'
             'Use -a "" to disable anti-cleavage filtering.',
    )
    parser.add_argument(
        '--cleavage_position',
        '-p',
        dest='cleavage_position',
        default=None,
        choices=['c', 'n'],
        help='Manual cleavage position. Overrides --enzyme.\n'
             'c = cut after the site, n = cut before the site.',
    )
    parser.add_argument(
        '--min_peptide_length',
        '-l',
        dest='min_len',
        default=6,
        type=int,
        help='Set minimum peptide length. Default = 6',
    )
    parser.add_argument(
        '--max_peptide_length',
        '-M',
        dest='max_len',
        default=40,
        type=int,
        help='Set maximum peptide length. Applied to all output modes. Default = 40',
    )
    parser.add_argument(
        '--missed-cleavages',
        '-L',
        dest='miss_cleavage',
        default=2,
        type=int,
        help='Maximum missed cleavages. Default = 2',
    )
    parser.add_argument(
        '--isobaric',
        default=False,
        action='store_true',
        help='Normalize I -> L before digestion.\n'
             'Default is to keep I and L distinct in this standalone CLI.',
    )
    parser.add_argument(
        '--initiator-methionine',
        choices=['both', 'remove', 'retain'],
        default='both',
        help='Handle an N-terminal M before digestion. Default = both.\n'
             'both = digest retained and removed forms; remove = digest without M;\n'
             'retain = digest the original sequence only.',
    )
    parser.add_argument(
        '--output-format',
        choices=['tsv', 'peptide', 'fasta'],
        default='tsv',
        help='Output format.\n'
             'tsv     -> protein_id\\tPEPTIDE\\tSTART\\tEND\n'
             'peptide -> unique PEPTIDE\n'
             'fasta   -> unique >protein_id_1\\nPEPTIDE',
    )
    parser.add_argument(
        '--header-template',
        default='{protein_id}_{index}',
        help='Header template for --output-format fasta.\n'
             'Placeholders: {protein_id}, {index}, {header}\n'
             'Examples:\n'
             '  {protein_id}_{index}\n'
             '  {protein_id}|pep{index}\n'
             '  pep{index}',
    )
    return parser


def _resolve_cli_cleavage(args):
    rule = resolve_cleavage_rule(
        enzyme=args.enzyme,
        cleavage_sites=args.cleavage_sites,
        anti_cleavage_sites=args.anti_cleavage_sites,
        cleavage_position=args.cleavage_position,
    )
    args.enzyme = rule['enzyme']
    args.cleavage_sites = rule['sites']
    args.anti_cleavage_sites = rule['no']
    args.cleavage_position = rule['pos']
    return args


def _iter_cli_results(args):
    if args.sequence:
        yield 'sequence', digest_sequence_occurrences(
            args.sequence,
            enzyme=args.enzyme,
            cleavage_sites=args.cleavage_sites,
            anti_cleavage_sites=args.anti_cleavage_sites,
            cleavage_position=args.cleavage_position,
            miss_cleavage=args.miss_cleavage,
            min_len=args.min_len,
            max_len=args.max_len,
            isobaric=args.isobaric,
            initiator_methionine=args.initiator_methionine,
        )
        return

    for file_fasta in args.fasta:
        yield from iter_digest_records(
            file_fasta,
            enzyme=args.enzyme,
            cleavage_sites=args.cleavage_sites,
            anti_cleavage_sites=args.anti_cleavage_sites,
            cleavage_position=args.cleavage_position,
            miss_cleavage=args.miss_cleavage,
            min_len=args.min_len,
            max_len=args.max_len,
            isobaric=args.isobaric,
            initiator_methionine=args.initiator_methionine,
        )


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if bool(args.sequence) == bool(args.fasta):
        parser.error('provide either FASTA input or --sequence')
    try:
        args = _resolve_cli_cleavage(args)
        write_digest_results(
            _iter_cli_results(args),
            sys.stdout,
            output_format=args.output_format,
            header_template=args.header_template,
        )
    except ValueError as exc:
        parser.error(str(exc))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
