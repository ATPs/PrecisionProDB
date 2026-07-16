import argparse
from bisect import bisect_left
import gzip
import re
from string import Formatter
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
HEADER_TEMPLATE_FIELDS = {'protein_id', 'index', 'header'}
DESCRIPTION = """Digest protein sequences from FASTA input or a direct sequence string.

Input modes:
  1. FASTA file(s):
       python src/precisionprodb/digestion.py proteins.fa proteins2.fa.gz
  2. One sequence:
       python src/precisionprodb/digestion.py --sequence AKRPQK

Cleavage rules:
  - Default behavior is the same as --enzyme Trypsin --missed-cleavages 2.
  - Use --enzyme for Comet-style presets.
  - Use --enzyme2 to add a second enzyme; cleavage positions are merged.
  - Use --cleavage-regex/--cleavage-regex2 for custom zero-width cleavage boundaries.
  - Use -c/-a/-p to override the preset manually.

Cleavage specificity:
  - full: both peptide termini must follow the cleavage rule (default).
  - n-specific: the N terminus follows the rule; the C terminus may be truncated.
  - c-specific: the C terminus follows the rule; the N terminus may be truncated.
  - semi: union of n-specific and c-specific results, with duplicates removed.

I/L handling in this standalone CLI:
  - Default: keep I and L distinct.
  - --isobaric: normalize I -> L before digestion.

Initiator methionine handling:
  - --initiator-methionine controls whether an N-terminal M is retained,
    removed, or both forms are digested (default: both).

Peptide attachments:
  - --nterm-attach and --cterm-attach add fixed strings after digestion.
  - Length filtering and TSV coordinates refer to the core digested peptide.

Stop-codon handling:
  - Internal "*" is always treated as a hard breakpoint.
  - Output peptides never include "*".
  - Internal "*" still counts in TSV peptide coordinates.

Output behavior:
  - TSV always includes protein_id, peptide, start, and end.
  - protein_id is the first whitespace-delimited token from the FASTA header.
  - FASTA records and peptide occurrences are streamed by the CLI.
  - peptide and fasta outputs are de-duplicated by peptide sequence.
"""
EPILOG = """Examples:
  Single-sequence digest:
    python src/precisionprodb/digestion.py --sequence AKRPQK -l 2

  FASTA digestion with missed cleavages:
    python src/precisionprodb/digestion.py proteins.fa --enzyme Trypsin --missed-cleavages 2 -l 6 -M 40

  Two-enzyme digestion:
    python src/precisionprodb/digestion.py proteins.fa --enzyme Trypsin --enzyme2 Lys_C

  Custom regex boundaries:
    python src/precisionprodb/digestion.py --sequence AKRPQK --cleavage-regex '(?<=[KR])(?!P)'
    python src/precisionprodb/digestion.py --sequence AAENLYFQG --cleavage-regex '(?<=ENLYFQ)(?=[GS])'

  Semi-specific digestion:
    python src/precisionprodb/digestion.py proteins.fa --specificity semi

  Gzipped FASTA input:
    python src/precisionprodb/digestion.py proteins.fa.gz --output-format peptide

  Remove the initiator methionine before digestion:
    python src/precisionprodb/digestion.py --sequence MAKRPQK --initiator-methionine remove

  Add fixed terminal sequences to every output peptide:
    python src/precisionprodb/digestion.py --sequence PEPTIDE --enzyme No_cut -l 1 --nterm-attach SY --cterm-attach GG

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

Comet-style enzyme presets for --enzyme NAME and --enzyme2 NAME:
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
  - No_cut with full specificity emits the full sequence as one peptide;
    partial specificity emits prefixes and/or suffixes.
  - Manual -c/-a/-p values override --enzyme only; --enzyme2 uses its preset.
  - With --enzyme2, cleavage positions from both enzymes are merged and duplicates count once.
  - Cleavage regexes must match zero-width amino-acid boundaries; they do not remove residues.
  - Internal "*" always splits the protein before digestion and is never kept in output peptides.
  - TSV coordinates are 1-based inclusive and internal "*" counts in the parent protein position.
  - peptide and fasta outputs keep only the first occurrence of each peptide sequence.
  - TSV output does not retain all generated peptides; peptide and fasta output retain
    a set of unique peptide sequences because those formats require de-duplication.
  - --initiator-methionine both emits retained and N-terminal-M-removed forms;
    shared peptide occurrences are emitted once.
  - --specificity semi is the de-duplicated union of n-specific and c-specific results.
  - Terminal attachments are added after digestion and do not affect cleavage,
    missed-cleavage counting, length filtering, or TSV coordinates.
  - This PrecisionProDB version is pure Python and does not use fast_digest/Cython.
"""

def open_text_file(filename):
    """Open a plain-text or gzipped text file from a string or path-like input."""
    filename = str(filename)
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


def _compile_cleavage_regex(cleavage_regex):
    """Compile a custom cleavage-boundary regex."""
    try:
        return re.compile(cleavage_regex)
    except (TypeError, re.error) as exc:
        raise ValueError('invalid cleavage regex: {}'.format(exc))


def resolve_cleavage_rule(enzyme=DEFAULT_ENZYME, cleavage_sites=None, anti_cleavage_sites=None, cleavage_position=None, cleavage_regex=None):
    """Resolve preset and manual enzyme settings into one cleavage rule."""
    if cleavage_regex is not None and any(
        value is not None
        for value in [cleavage_sites, anti_cleavage_sites, cleavage_position]
    ):
        raise ValueError('cleavage_regex cannot be combined with cleavage_sites, anti_cleavage_sites, or cleavage_position')
    enzyme_name = resolve_enzyme_name(enzyme) or DEFAULT_ENZYME
    preset = ENZYME_PRESETS[enzyme_name]
    sites = preset['sites'] if cleavage_sites is None else cleavage_sites
    anti_sites = preset['no'] if anti_cleavage_sites is None else anti_cleavage_sites
    position = preset['pos'] if cleavage_position is None else cleavage_position
    if position not in ['c', 'n']:
        raise ValueError("cleavage_position must be 'c' or 'n'")
    rule = {
        'enzyme': enzyme_name,
        'sites': sites,
        'no': anti_sites,
        'pos': position,
    }
    if cleavage_regex is not None:
        rule['regex'] = _compile_cleavage_regex(cleavage_regex)
    return rule


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


def _regex_cleavage_boundaries(protein, cleavage_regex):
    """Return boundaries from zero-width regex matches plus protein termini."""
    boundaries = [0]
    for match in cleavage_regex.finditer(protein):
        if match.start() != match.end():
            raise ValueError(
                'cleavage regex must use zero-width matches; matched residues at {}-{}'.format(
                    match.start(),
                    match.end(),
                )
            )
        if match.start() != boundaries[-1]:
            boundaries.append(match.start())
    if boundaries[-1] != len(protein):
        boundaries.append(len(protein))
    return boundaries


def _rule_cleavage_boundaries(protein, rule):
    """Return boundaries for a preset/manual rule or a custom regex rule."""
    cleavage_regex = rule.get('regex')
    if cleavage_regex is not None:
        return _regex_cleavage_boundaries(protein, cleavage_regex)
    return _cleavage_boundaries(
        protein,
        sites=rule['sites'],
        pos=rule['pos'],
        no=rule['no'],
    )


def _is_cut_everywhere_resolved_rule(rule):
    """Return True only when a rule uses the Cut_everywhere preset behavior."""
    return (
        rule.get('regex') is None
        and _is_cut_everywhere_rule(rule['sites'], rule['pos'], rule['no'])
    )


def _merged_cleavage_boundaries(protein, rules):
    """Return the sorted union of cleavage boundaries from one or two rules."""
    boundaries = _rule_cleavage_boundaries(protein, rules[0])
    for rule in rules[1:]:
        additional_boundaries = _rule_cleavage_boundaries(protein, rule)
        merged = []
        first_index = 0
        second_index = 0
        while first_index < len(boundaries) or second_index < len(additional_boundaries):
            if second_index >= len(additional_boundaries):
                position = boundaries[first_index]
                first_index += 1
            elif first_index >= len(boundaries):
                position = additional_boundaries[second_index]
                second_index += 1
            elif boundaries[first_index] < additional_boundaries[second_index]:
                position = boundaries[first_index]
                first_index += 1
            elif additional_boundaries[second_index] < boundaries[first_index]:
                position = additional_boundaries[second_index]
                second_index += 1
            else:
                position = boundaries[first_index]
                first_index += 1
                second_index += 1
            if not merged or merged[-1] != position:
                merged.append(position)
        boundaries = merged
    return boundaries


def _segment_digest_occurrences(protein, sites='KR', pos='c', no='P', min_len=0, max_len=None, boundaries=None):
    """Yield peptide occurrences as ``(peptide, start_0based, end_0based_exclusive)``."""
    if boundaries is None:
        boundaries = _cleavage_boundaries(protein, sites=sites, pos=pos, no=no)
    for start, end in zip(boundaries, boundaries[1:]):
        peptide = protein[start:end]
        peptide_length = len(peptide)
        if not peptide or peptide_length < min_len:
            continue
        if max_len is not None and peptide_length > max_len:
            continue
        yield peptide, start, end


def _segment_cut_everywhere_occurrences(protein, min_len=0, max_len=None):
    """Yield all contiguous peptide occurrences for the requested length range."""
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
            yield protein[start:end], start, end


def _segment_digest_with_missed_cleavages_occurrences(protein, sites='KR', pos='c', no='P', miss_cleavage=2, peplen_min=6, peplen_max=40, boundaries=None):
    """Yield missed-cleavage occurrences as ``(peptide, start, end)``."""
    if boundaries is None:
        boundaries = _cleavage_boundaries(protein, sites=sites, pos=pos, no=no)
    max_missed = max(0, int(miss_cleavage))
    for missed_count in range(max_missed + 1):
        for start_index in range(len(boundaries) - missed_count - 1):
            start = boundaries[start_index]
            end = boundaries[start_index + missed_count + 1]
            peptide = protein[start:end]
            peptide_length = len(peptide)
            if peplen_min <= peptide_length <= peplen_max:
                yield peptide, start, end


def _segment_n_specific_occurrences(protein, sites='KR', pos='c', no='P', miss_cleavage=2, peplen_min=6, peplen_max=40, boundaries=None):
    """Yield occurrences with a specific N terminus and arbitrary C terminus."""
    if boundaries is None:
        boundaries = _cleavage_boundaries(protein, sites=sites, pos=pos, no=no)
    max_missed = max(0, int(miss_cleavage))
    min_length = max(1, peplen_min)
    for start_index in range(len(boundaries) - 1):
        start = boundaries[start_index]
        last_index = min(start_index + max_missed + 1, len(boundaries) - 1)
        min_end = start + min_length
        max_end = min(boundaries[last_index], start + peplen_max)
        for end in range(min_end, max_end + 1):
            yield protein[start:end], start, end


def _segment_c_specific_occurrences(protein, sites='KR', pos='c', no='P', miss_cleavage=2, peplen_min=6, peplen_max=40, boundaries=None):
    """Yield occurrences with an arbitrary N terminus and specific C terminus."""
    if boundaries is None:
        boundaries = _cleavage_boundaries(protein, sites=sites, pos=pos, no=no)
    max_missed = max(0, int(miss_cleavage))
    min_length = max(1, peplen_min)
    for end_index in range(1, len(boundaries)):
        end = boundaries[end_index]
        first_index = max(0, end_index - max_missed - 1)
        min_start = max(boundaries[first_index], end - peplen_max)
        max_start = end - min_length
        for start in range(min_start, max_start + 1):
            yield protein[start:end], start, end


def _segment_digest_with_specificity_occurrences(protein, sites='KR', pos='c', no='P', miss_cleavage=2, peplen_min=6, peplen_max=40, specificity='full', boundaries=None):
    """Yield peptide occurrences for the requested cleavage specificity."""
    if specificity == 'full':
        yield from _segment_digest_with_missed_cleavages_occurrences(
            protein,
            sites=sites,
            pos=pos,
            no=no,
            miss_cleavage=miss_cleavage,
            peplen_min=peplen_min,
            peplen_max=peplen_max,
            boundaries=boundaries,
        )
        return
    if specificity == 'n-specific':
        yield from _segment_n_specific_occurrences(
            protein,
            sites=sites,
            pos=pos,
            no=no,
            miss_cleavage=miss_cleavage,
            peplen_min=peplen_min,
            peplen_max=peplen_max,
            boundaries=boundaries,
        )
        return
    if specificity == 'c-specific':
        yield from _segment_c_specific_occurrences(
            protein,
            sites=sites,
            pos=pos,
            no=no,
            miss_cleavage=miss_cleavage,
            peplen_min=peplen_min,
            peplen_max=peplen_max,
            boundaries=boundaries,
        )
        return
    if specificity == 'semi':
        if boundaries is None:
            boundaries = _cleavage_boundaries(protein, sites=sites, pos=pos, no=no)
        for occurrence in _segment_n_specific_occurrences(
            protein,
            sites=sites,
            pos=pos,
            no=no,
            miss_cleavage=miss_cleavage,
            peplen_min=peplen_min,
            peplen_max=peplen_max,
            boundaries=boundaries,
        ):
            yield occurrence
        for occurrence in _segment_c_specific_occurrences(
            protein,
            sites=sites,
            pos=pos,
            no=no,
            miss_cleavage=miss_cleavage,
            peplen_min=peplen_min,
            peplen_max=peplen_max,
            boundaries=boundaries,
        ):
            boundary_index = bisect_left(boundaries, occurrence[1])
            if boundary_index >= len(boundaries) or boundaries[boundary_index] != occurrence[1]:
                yield occurrence
        return
    raise ValueError("specificity must be 'full', 'n-specific', 'c-specific', or 'semi'")


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


def iter_digest_sequence_occurrences(sequence, enzyme=DEFAULT_ENZYME, cleavage_sites=None, anti_cleavage_sites=None, cleavage_position=None, miss_cleavage=2, min_len=6, max_len=40, isobaric=False, initiator_methionine='both', specificity='full', enzyme2=None, nterm_attach='', cterm_attach='', cleavage_regex=None, cleavage_regex2=None):
    """Yield peptide occurrences from one protein sequence."""
    rule = resolve_cleavage_rule(
        enzyme=enzyme,
        cleavage_sites=cleavage_sites,
        anti_cleavage_sites=anti_cleavage_sites,
        cleavage_position=cleavage_position,
        cleavage_regex=cleavage_regex,
    )
    rules = [rule]
    if enzyme2 or cleavage_regex2 is not None:
        rules.append(
            resolve_cleavage_rule(
                enzyme=enzyme2 or DEFAULT_ENZYME,
                cleavage_regex=cleavage_regex2,
            )
        )
    if specificity not in ['full', 'n-specific', 'c-specific', 'semi']:
        raise ValueError("specificity must be 'full', 'n-specific', 'c-specific', or 'semi'")
    nterm_attach = str(nterm_attach or '')
    cterm_attach = str(cterm_attach or '')
    normalized_sequence = normalize_sequence(sequence, isobaric=isobaric)
    cut_everywhere = any(_is_cut_everywhere_resolved_rule(item) for item in rules)
    variant_mode = initiator_methionine
    if cut_everywhere and variant_mode == 'both' and normalized_sequence.startswith('M'):
        variant_mode = 'retain'
    deduplicate_variants = (
        variant_mode == 'both'
        and normalized_sequence.startswith('M')
    )
    seen_positions = set() if deduplicate_variants else None
    position_multiplier = len(normalized_sequence) + 1
    for sequence_offset, sequence_variant in iter_initiator_methionine_variants(
        normalized_sequence,
        initiator_methionine=variant_mode,
    ):
        for segment_offset, segment in iter_stop_segments(sequence_variant):
            if cut_everywhere:
                segment_occurrences = _segment_cut_everywhere_occurrences(
                    segment,
                    min_len=min_len,
                    max_len=max_len,
                )
            else:
                boundaries = _merged_cleavage_boundaries(segment, rules)
                segment_occurrences = _segment_digest_with_specificity_occurrences(
                    segment,
                    sites=rule['sites'],
                    pos=rule['pos'],
                    no=rule['no'],
                    miss_cleavage=miss_cleavage,
                    peplen_min=min_len,
                    peplen_max=max_len,
                    specificity=specificity,
                    boundaries=boundaries,
                )
            for peptide, start, end in segment_occurrences:
                occurrence = (
                    nterm_attach + peptide + cterm_attach,
                    sequence_offset + segment_offset + start + 1,
                    sequence_offset + segment_offset + end,
                )
                if seen_positions is not None:
                    position_key = occurrence[1] * position_multiplier + occurrence[2]
                    if position_key in seen_positions:
                        continue
                    seen_positions.add(position_key)
                yield occurrence


def digest_sequence_occurrences(sequence, enzyme=DEFAULT_ENZYME, cleavage_sites=None, anti_cleavage_sites=None, cleavage_position=None, miss_cleavage=2, min_len=6, max_len=40, isobaric=False, initiator_methionine='both', specificity='full', enzyme2=None, nterm_attach='', cterm_attach='', cleavage_regex=None, cleavage_regex2=None):
    """Digest one protein sequence and return a list of peptide occurrences."""
    return list(
        iter_digest_sequence_occurrences(
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
            specificity=specificity,
            enzyme2=enzyme2,
            nterm_attach=nterm_attach,
            cterm_attach=cterm_attach,
            cleavage_regex=cleavage_regex,
            cleavage_regex2=cleavage_regex2,
        )
    )


def digest_sequence(sequence, enzyme=DEFAULT_ENZYME, cleavage_sites=None, anti_cleavage_sites=None, cleavage_position=None, miss_cleavage=2, min_len=6, max_len=40, isobaric=False, initiator_methionine='both', specificity='full', enzyme2=None, nterm_attach='', cterm_attach='', cleavage_regex=None, cleavage_regex2=None):
    """Digest one protein sequence with missed cleavages."""
    return [
        peptide
        for peptide, _start, _end in iter_digest_sequence_occurrences(
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
            specificity=specificity,
            enzyme2=enzyme2,
            nterm_attach=nterm_attach,
            cterm_attach=cterm_attach,
            cleavage_regex=cleavage_regex,
            cleavage_regex2=cleavage_regex2,
        )
    ]


def iter_digest_records(file_path, enzyme=DEFAULT_ENZYME, cleavage_sites=None, anti_cleavage_sites=None, cleavage_position=None, miss_cleavage=2, min_len=6, max_len=40, isobaric=False, initiator_methionine='both', specificity='full', enzyme2=None, nterm_attach='', cterm_attach='', cleavage_regex=None, cleavage_regex2=None):
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
            specificity=specificity,
            enzyme2=enzyme2,
            nterm_attach=nterm_attach,
            cterm_attach=cterm_attach,
            cleavage_regex=cleavage_regex,
            cleavage_regex2=cleavage_regex2,
        )


def iter_digest_records_streaming(file_path, enzyme=DEFAULT_ENZYME, cleavage_sites=None, anti_cleavage_sites=None, cleavage_position=None, miss_cleavage=2, min_len=6, max_len=40, isobaric=False, initiator_methionine='both', specificity='full', enzyme2=None, nterm_attach='', cterm_attach='', cleavage_regex=None, cleavage_regex2=None):
    """Yield FASTA records with streaming peptide-occurrence iterators."""
    for header, sequence in read_fasta_records(file_path):
        yield header, iter_digest_sequence_occurrences(
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
            specificity=specificity,
            enzyme2=enzyme2,
            nterm_attach=nterm_attach,
            cterm_attach=cterm_attach,
            cleavage_regex=cleavage_regex,
            cleavage_regex2=cleavage_regex2,
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
        for _literal, field_name, _format_spec, _conversion in Formatter().parse(template):
            if field_name is not None and field_name not in HEADER_TEMPLATE_FIELDS:
                raise ValueError('unsupported header template field: {}'.format(field_name))
        return template.format(protein_id=protein_id, index=index, header=header_text)
    except (AttributeError, IndexError, KeyError, ValueError) as exc:
        raise ValueError('invalid header template: {}'.format(exc))


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
        '--enzyme2',
        default='',
        help='Optional second Comet-style enzyme preset.\n'
             'Its cleavage positions are merged with --enzyme.',
    )
    parser.add_argument(
        '--cleavage-regex',
        default=None,
        help='Custom zero-width regex for first-rule cleavage boundaries.\n'
             'Overrides --enzyme and cannot be combined with -c/-a/-p.',
    )
    parser.add_argument(
        '--cleavage-regex2',
        default=None,
        help='Custom zero-width regex for second-rule cleavage boundaries.\n'
             'Overrides --enzyme2 and may be used without --enzyme2.',
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
        '--specificity',
        choices=['full', 'n-specific', 'c-specific', 'semi'],
        default='full',
        help='Cleavage specificity. Default = full.\n'
             'full = both termini specific; n-specific = arbitrary C terminus;\n'
             'c-specific = arbitrary N terminus; semi = union of both partial modes.',
    )
    parser.add_argument(
        '--nterm-attach',
        default='',
        help='Fixed string to prepend to every digested peptide. Default = empty.',
    )
    parser.add_argument(
        '--cterm-attach',
        default='',
        help='Fixed string to append to every digested peptide. Default = empty.',
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
    if args.min_len < 1:
        raise ValueError('min_peptide_length must be at least 1')
    if args.max_len < 1:
        raise ValueError('max_peptide_length must be at least 1')
    if args.min_len > args.max_len:
        raise ValueError('min_peptide_length cannot exceed max_peptide_length')
    if args.miss_cleavage < 0:
        raise ValueError('missed-cleavages must be at least 0')
    resolve_cleavage_rule(
        enzyme=args.enzyme,
        cleavage_sites=args.cleavage_sites,
        anti_cleavage_sites=args.anti_cleavage_sites,
        cleavage_position=args.cleavage_position,
        cleavage_regex=args.cleavage_regex,
    )
    if args.enzyme2 or args.cleavage_regex2 is not None:
        resolve_cleavage_rule(
            enzyme=args.enzyme2 or DEFAULT_ENZYME,
            cleavage_regex=args.cleavage_regex2,
        )
    return args


def _iter_cli_results(args):
    if args.sequence:
        yield 'sequence', iter_digest_sequence_occurrences(
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
            specificity=args.specificity,
            enzyme2=args.enzyme2,
            nterm_attach=args.nterm_attach,
            cterm_attach=args.cterm_attach,
            cleavage_regex=args.cleavage_regex,
            cleavage_regex2=args.cleavage_regex2,
        )
        return

    for file_fasta in args.fasta:
        yield from iter_digest_records_streaming(
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
            specificity=args.specificity,
            enzyme2=args.enzyme2,
            nterm_attach=args.nterm_attach,
            cterm_attach=args.cterm_attach,
            cleavage_regex=args.cleavage_regex,
            cleavage_regex2=args.cleavage_regex2,
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
