import gzip
import io
import glob
import os
from multiprocessing import Pool
import numpy as np
import shutil
import subprocess

if __package__:
    from .args import add_argument_set
    from .runstate import StageCache
else:
    from args import add_argument_set
    from runstate import StageCache

try:
    import tqdm
except:
    print('Cannot import tqdm. Will not use tqdm.')
    tqdm = False


def _import_pandas():
    import pandas as pd
    return pd


CHROMOSOMES = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'}

def openFile(filename):
    '''open text or gzipped text file
    '''
    if filename.endswith('.gz'):
        return gzip.open(filename,'rt')
    return open(filename)


def strip_vcf_suffix(filename):
    '''Return a basename suitable for sample naming by removing common VCF suffixes.'''
    base_name = os.path.basename(filename)
    for suffix in ('.vcf.gz', '.vcf', '.gz'):
        if base_name.lower().endswith(suffix):
            return base_name[:-len(suffix)]
    return os.path.splitext(base_name)[0]


def is_manifest_file(file_input):
    '''Return True when an existing TSV-like file uses filepath as the first header column.'''
    if not os.path.exists(file_input):
        return False
    try:
        with openFile(file_input) as fo:
            header = fo.readline().rstrip('\n\r').split('\t')
    except Exception:
        return False
    return len(header) > 0 and header[0] == 'filepath'


def get_vcf_samples(file_vcf):
    '''Read the #CHROM header and return sample names from a VCF/VCF.GZ file.'''
    with openFile(file_vcf) as fo:
        for line in fo:
            if line.startswith('#CHROM'):
                columns = line.rstrip('\n\r').split('\t')
                return columns[9:]
    raise ValueError(f'#CHROM header line not found in {file_vcf}')


def read_vcf_manifest(file_manifest):
    '''Parse a VCF manifest and resolve one unique output sample name per row.'''
    pd = _import_pandas()
    df_manifest = pd.read_csv(file_manifest, sep='\t', keep_default_na=False)
    if df_manifest.shape[1] == 0 or df_manifest.columns[0] != 'filepath':
        raise ValueError('manifest TSV must have header and first column named filepath')
    if df_manifest.shape[0] == 0:
        raise ValueError('manifest TSV has no input rows')

    rows = []
    for row_index, row in df_manifest.iterrows():
        file_vcf = str(row['filepath']).strip()
        if file_vcf == '':
            raise ValueError(f'empty filepath in manifest row {row_index + 2}')
        if not os.path.exists(file_vcf):
            raise FileNotFoundError(f'manifest filepath not found: {file_vcf}')
        if not (file_vcf.lower().endswith('.vcf') or file_vcf.lower().endswith('.vcf.gz')):
            raise ValueError(f'manifest filepath must point to a VCF/VCF.GZ file in this version: {file_vcf}')
        samples = get_vcf_samples(file_vcf)
        if not samples:
            raise ValueError(f'no sample columns found in {file_vcf}')
        sample = str(row['sample']).strip() if 'sample' in df_manifest.columns else ''
        if sample == '':
            sample = samples[0]
            print(f'Warning: sample not provided for {file_vcf}; using first sample {sample}')
        if sample not in samples:
            raise ValueError(f'sample {sample} not found in {file_vcf}')
        sample_col = 9 + samples.index(sample)
        name_use = str(row['name_use']).strip() if 'name_use' in df_manifest.columns else ''
        if name_use == '':
            name_use = sample
        rows.append({'filepath': file_vcf, 'sample': sample, 'sample_col': sample_col, 'name_use': name_use})

    names = [row['name_use'] for row in rows]
    duplicate_names = {name for name in names if names.count(name) > 1}
    if duplicate_names:
        print('Warning: duplicated sample names in manifest output:', ','.join(sorted(duplicate_names)))
        print('Warning: using VCF file names as sample names for duplicated entries.')
        for row in rows:
            if row['name_use'] in duplicate_names:
                row['name_use'] = strip_vcf_suffix(row['filepath'])

    names = [row['name_use'] for row in rows]
    duplicate_names = {name for name in names if names.count(name) > 1}
    if duplicate_names:
        raise ValueError('duplicated sample names remain after using file names; please provide unique name_use values: ' + ','.join(sorted(duplicate_names)))

    return rows


def openVCFStreamFast(filename, threads=1):
    '''Open VCF text stream, using external parallel gzip decompression when available.'''
    if not filename.endswith('.gz'):
        return open(filename), None

    decompress_threads = max(1, min(4, int(threads) if threads else 1))
    if shutil.which('pigz'):
        cmd = ['pigz', '-dc', '-p', str(decompress_threads), filename]
    elif shutil.which('bgzip'):
        cmd = ['bgzip', '-@', str(decompress_threads), '-dc', filename]
    else:
        return gzip.open(filename, 'rt'), None

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return io.TextIOWrapper(proc.stdout), proc


def closeVCFStreamFast(fo, proc=None):
    '''Close a VCF stream opened by openVCFStreamFast and check decompressor status.'''
    fo.close()
    if proc is not None:
        _, stderr = proc.communicate()
        if proc.returncode not in (0, None):
            raise RuntimeError(stderr.decode(errors='replace').strip())


def iter_vcf_records(file_vcf, threads=1):
    '''Yield non-header VCF records from a text or gzipped VCF.'''
    fo, proc = openVCFStreamFast(file_vcf, threads=threads)
    try:
        for line in fo:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                break
        for line in fo:
            yield line
    finally:
        closeVCFStreamFast(fo, proc)


def info_field_value(info_str, field):
    '''Return a single INFO field value without parsing the whole INFO column.'''
    start = info_str.find(field + '=')
    if start == -1:
        return None
    if start > 0 and info_str[start - 1] != ';':
        start = info_str.find(';' + field + '=')
        if start == -1:
            return None
        start += 1
    value_start = start + len(field) + 1
    value_end = info_str.find(';', value_start)
    if value_end == -1:
        value_end = len(info_str)
    return info_str[value_start:value_end]


def info_field_passes_threshold(info_str, field, threshold):
    '''Return True if any comma-separated INFO value reaches the threshold.'''
    value = info_field_value(info_str, field)
    if value is None:
        return None
    for one_value in value.split(','):
        try:
            if float(one_value) >= threshold:
                return True
        except ValueError:
            continue
    return False


def _resolve_sample_selection(columns, individual_input, file_vcf):
    """Return selected sample names and VCF column indices.

    Explicit sample requests are all-or-nothing: a typo must not silently turn
    a sample-specific conversion into an all-variants conversion.
    """
    all_samples = columns[9:] if len(columns) > 9 else []
    if individual_input == 'ALL_VARIANTS':
        individuals = []
    elif individual_input == 'ALL_SAMPLES':
        individuals = all_samples
    elif individual_input is None:
        if all_samples:
            individuals = [all_samples[0]]
            print(f"No sample specified, using first sample: {individuals[0]}")
        else:
            individuals = []
            print(f"No sample columns found in {file_vcf}; processing all variants.")
    else:
        if isinstance(individual_input, str):
            requested = [sample.strip() for sample in individual_input.split(',')]
        else:
            requested = [str(sample).strip() for sample in individual_input]
        requested = list(dict.fromkeys(sample for sample in requested if sample))
        if not requested:
            raise ValueError(f"No sample names were requested for {file_vcf}")
        missing = [sample for sample in requested if sample not in all_samples]
        if missing:
            raise ValueError(
                f"Requested sample(s) not found in {file_vcf}: {', '.join(missing)}"
            )
        individuals = requested

    individual_col = [9 + all_samples.index(sample) for sample in individuals]
    return individuals, individual_col


def _decode_genotype(sample_value, format_field, n_alternatives, context):
    """Decode one VCF sample field as two allele indices.

    Missing alleles are treated as reference, and haploid calls receive a
    reference allele in the second slot to preserve the mutation TSV schema.
    """
    format_keys = format_field.split(':') if format_field else []
    if 'GT' not in format_keys:
        raise ValueError(f"GT field is missing for {context} (FORMAT={format_field!r})")
    gt_index = format_keys.index('GT')
    sample_fields = sample_value.split(':')
    if gt_index >= len(sample_fields):
        raise ValueError(f"GT value is missing for {context}")
    genotype = sample_fields[gt_index]
    if genotype == '':
        raise ValueError(f"Empty GT value for {context}")

    if '|' in genotype and '/' in genotype:
        raise ValueError(f"Mixed phased and unphased separators in GT {genotype!r} for {context}")
    if '|' in genotype:
        alleles = genotype.split('|')
    elif '/' in genotype:
        alleles = genotype.split('/')
    else:
        alleles = [genotype]

    if len(alleles) > 2:
        raise ValueError(f"More than two alleles in GT {genotype!r} for {context}")
    if len(alleles) == 1:
        alleles.append('0')

    decoded = []
    for allele in alleles:
        if allele == '.':
            decoded.append(0)
            continue
        if not allele.isdigit():
            raise ValueError(f"Invalid allele {allele!r} in GT {genotype!r} for {context}")
        allele_index = int(allele)
        if allele_index > n_alternatives:
            raise ValueError(
                f"GT allele index {allele_index} exceeds {n_alternatives} ALT allele(s) for {context}"
            )
        decoded.append(allele_index)
    return decoded


def _decode_record_genotypes(es, individual_col, file_vcf=None, individual_names=None):
    """Decode selected genotypes from split VCF fields using one shared rule."""
    chromosome = es[0] if es else '<unknown chromosome>'
    position = es[1] if len(es) > 1 else '<unknown position>'
    context_base = f"{file_vcf or '<VCF>'}:{chromosome}:{position}"
    alternatives = es[4].split(',') if len(es) > 4 else []
    format_field = es[8] if len(es) > 8 else ''
    genotypes = []
    for selected_index, column_index in enumerate(individual_col):
        if column_index >= len(es):
            sample_label = (
                individual_names[selected_index]
                if individual_names and selected_index < len(individual_names)
                else f'column {column_index + 1}'
            )
            raise ValueError(f"Missing sample field for {sample_label} at {context_base}")
        sample_label = (
            individual_names[selected_index]
            if individual_names and selected_index < len(individual_names)
            else f'column {column_index + 1}'
        )
        context = f"sample {sample_label} at {context_base}"
        genotypes.append(
            _decode_genotype(es[column_index], format_field, len(alternatives), context)
        )
    return genotypes

def getMutationsFromVCF(file_vcf, outprefix = None, individual=None, filter_PASS = True, chromosome_only = True):
    '''return two dataframe of mutations.
    if outprefix is not None, write the two dataframe to outprfex +'_1/2.tsv' two files as mutation file
    file_vcf is the vcf file, it could be a gzip file
    if individual is None, work on the first individual in the vcf file
    default, only work on variants in chromosome
    '''
    pd = _import_pandas()
    if ',' in file_vcf:
        files_vcf = file_vcf.split(',')
    elif os.path.exists(file_vcf):
        files_vcf = [file_vcf]
    else:
        files_vcf = glob.glob(file_vcf)
        if len(files_vcf) == 0:
            print(file_vcf, 'not found.')
            return None
    
    ls_keep = ['chr\tpos\tref\talt1\talt2\n']

    for file_vcf in files_vcf:
        with openFile(file_vcf) as fo:
            header_line = None
            for line in fo:
                if line.startswith('##'):
                    continue
                header_line = line
                break
            if header_line is None or not header_line.startswith('#CHROM'):
                raise ValueError(f"#CHROM header line not found in {file_vcf}")

            columns = header_line.rstrip('\n\r').split('\t')
            selected, individual_col = _resolve_sample_selection(columns, individual, file_vcf)
            if len(selected) != 1:
                raise ValueError(
                    "legacy VCF conversion requires exactly one sample; use "
                    "convertVCF2MutationComplex for multiple samples or ALL_VARIANTS"
                )
            print('using sample', selected[0], 'from', file_vcf)

            for line in fo:
                if not line.strip():
                    continue
                es = line.rstrip('\n\r').split('\t')
                if len(es) < 8:
                    raise ValueError(f"Malformed VCF record with fewer than 8 columns in {file_vcf}: {line.rstrip()}")
                chromosome, position, reference, alternatives, filter_value = es[0], es[1], es[3], es[4], es[6]

                if chromosome_only and chromosome not in CHROMOSOMES:
                    continue
                if filter_PASS and filter_value != 'PASS':
                    continue

                genotype = _decode_record_genotypes(es, individual_col, file_vcf, selected)[0]
                if all(allele == 0 for allele in genotype):
                    continue

                alleles = [reference] + alternatives.split(',')
                alternative1, alternative2 = alleles[genotype[0]], alleles[genotype[1]]
                ls_keep.append(
                    f'{chromosome}\t{position}\t{reference}\t{alternative1}\t{alternative2}\n'
                )
    
    df = pd.read_csv(io.StringIO(''.join(ls_keep)), sep='\t',low_memory=False)
    print('before QC, sites with mutations:', df.shape[0],'of which, homozyous sites:', df[df['alt1'] == df['alt2']].shape[0])

    df1 = df[['chr', 'pos', 'ref', 'alt1']].copy()
    df2 = df[['chr', 'pos', 'ref', 'alt2']].copy()
    column_names = ['chr', 'pos', 'ref', 'alt']
    df1.columns = column_names
    df2.columns = column_names

    # filter dataframe
    df1 = df1[(df1['ref'] != df1['alt']) & (df1['alt'] != '*')]
    df2 = df2[(df2['ref'] != df2['alt']) & (df2['alt'] != '*')]
    df1 = df1.drop_duplicates(subset=['chr','pos','ref'], keep='first')
    df2 = df2.drop_duplicates(subset=['chr','pos','ref'], keep='first')
    print('after QC, number of variants in two dataframes are:', df1.shape[0], df2.shape[0])

    # save results
    if outprefix is not None:
        fout1 = outprefix + '_1.tsv'
        fout2 = outprefix + '_2.tsv'
        df1.to_csv(fout1, sep='\t', index=None)
        df2.to_csv(fout2, sep='\t', index=None)
        return [fout1, fout2]

    return df1, df2


def split_complex_variant(ls_to_write):
    '''
    ls_to_write = [chromosome, position, reference, alternative, '\t'.join(str(e) for e in ls_alternatives)]
    or ls_to_write = [chromosome, position, reference, alternative]
    deal with complex variants. usually it should be like chr1-111-A-T, chr1-111-AT-A, chr1-111-A-AT. but some times, it could be like chr1-111-GG-TT
    
    Args:
        ls_to_write (list): A list representing a variant row. Expected format:
            [chromosome, position, reference, alternative, optional_col1, optional_col2, ...]
            Example: ['chr1', '100', 'A', 'T']
            Example: ['chr1', '200', 'G', 'TT', '0', '1']
            Example: ['chr1', '300', 'GG', 'TT', '1', '1']

    Returns:
        list[list]: A list containing one or more lists, where each inner list
                    represents a decomposed variant row ready for output.
                    Format: [chromosome, position, reference, alternative, optional_col1, ...]

    '''
    chromosome, position, reference, alternative = ls_to_write[:4]
    
    len_ref = len(reference)
    len_alt = len(alternative)
    # If it's a simple SNP, insertion, or deletion (at least one allele has length 1)
    # No decomposition is needed according to the original logic's effective outcome.
    if len_ref == 1 or len_alt == 1:
        # Return the original variant wrapped in a list
        return [ls_to_write]
    
    
    position = int(position)
    ls = []
    n = min(len(reference), len(alternative))
    for i in range(n-1):
        ref = reference[i]
        alt = alternative[i]
        if ref != alt:
            ls.append([chromosome, position + i, ref, alt] + ls_to_write[4:])
    ref = reference[n-1:]
    alt = alternative[n-1:]
    ls.append([chromosome, position + n-1, ref, alt] + ls_to_write[4:])
    return ls

def get_writing_string_for_complex_variant(ls_to_write):
    '''
    ls_to_write = [chromosome, position, reference, alternative, '\t'.join(str(e) for e in ls_alternatives)]
    or ls_to_write = [chromosome, position, reference, alternative]
    deal with complex variants. usually it should be like chr1-111-A-T, chr1-111-AT-A, chr1-111-A-AT. but some times, it could be like chr1-111-GG-TT
    
    '''
    split_variants = split_complex_variant(ls_to_write)
    ls_txt = ['\t'.join(map(str, variant_row)) + '\n' for variant_row in split_variants]
    return ''.join(ls_txt)


def processOneLineOfVCF(line, 
                        individual_col, 
                        chromosome_only=True, 
                        filter_PASS=True, 
                        info_field=None, 
                        info_field_thres=None, 
                        chromosomes=CHROMOSOMES,
                        file_vcf=None,
                        individual_names=None,
                        ):
    '''Process one line of VCF file and return formatted mutation information
    
    Args:
        line (str): One line from VCF file
        individual_col (list): List of column indices for individuals/samples
        chromosome_only (bool): If True, only process chromosomes (skip scaffolds)
        filter_PASS (bool): If True, only process variants with FILTER=PASS
        info_field (str): INFO field key to filter by (e.g. 'AF')
        info_field_thres (float): Threshold for info_field filtering
        chromosomes (set): Set of valid chromosome names
        file_vcf (str): VCF file path for error reporting
        
    Returns:
        str: Formatted mutation line(s) or empty string if variant should be skipped
    '''
    es = line.rstrip('\n\r').split('\t')
    if len(es) < 8:
        raise ValueError(f'Malformed VCF line with fewer than 8 columns in {file_vcf}: {line.rstrip()}')
    chromosome, position, reference, alternatives, FILTER = es[0], es[1], es[3], es[4], es[6]

    if chromosome_only:
        # skip if not chromosome but some scaffolds
        if chromosome not in chromosomes:
            return ''
    
    if filter_PASS:
        # skip if not "PASS"
        if FILTER != "PASS":
            return ''

    # filter by info field
    if info_field:
        if info_field_thres:
            info_str = es[7]
            passes_threshold = info_field_passes_threshold(info_str, info_field, info_field_thres)
            if passes_threshold is None:
                print('waring! INFO field', info_field, 'not found in', file_vcf, line)
            elif not passes_threshold:
                return ''
    alternatives = alternatives.split(',')
    GTs = _decode_record_genotypes(es, individual_col, file_vcf, individual_names)
    if GTs and all(all(allele == 0 for allele in genotype) for genotype in GTs):
        return ''
    
    alleles = [reference] + alternatives
    ls_new_line = []
    for alternative_index in range(len(alternatives)):
        alternative = alternatives[alternative_index]
        if alternative == '*':
            continue
        ls_alternatives = []
        for GT in GTs:
            for i in GT:
                if i == alternative_index + 1:
                    ls_alternatives.append(1)
                else:
                    ls_alternatives.append(0)
        if ls_alternatives:
            ls_to_write = [chromosome, position, reference, alternative, '\t'.join(str(e) for e in ls_alternatives)]
        else:
            ls_to_write = [chromosome, position, reference, alternative]
        new_line = get_writing_string_for_complex_variant(ls_to_write)
        ls_new_line.append(new_line)
    
    return ''.join(ls_new_line)

_WORKER_INDIVIDUAL_COL = []
_WORKER_CHROMOSOME_ONLY = True
_WORKER_FILTER_PASS = True
_WORKER_INFO_FIELD = None
_WORKER_INFO_FIELD_THRES = None
_WORKER_FILE_VCF = None
_WORKER_INDIVIDUAL_NAMES = None


def init_vcf_worker(individual_col, chromosome_only, filter_PASS, info_field, info_field_thres, file_vcf, individual_names=None):
    global _WORKER_INDIVIDUAL_COL
    global _WORKER_CHROMOSOME_ONLY
    global _WORKER_FILTER_PASS
    global _WORKER_INFO_FIELD
    global _WORKER_INFO_FIELD_THRES
    global _WORKER_FILE_VCF
    global _WORKER_INDIVIDUAL_NAMES
    _WORKER_INDIVIDUAL_COL = individual_col
    _WORKER_CHROMOSOME_ONLY = chromosome_only
    _WORKER_FILTER_PASS = filter_PASS
    _WORKER_INFO_FIELD = info_field
    _WORKER_INFO_FIELD_THRES = info_field_thres
    _WORKER_FILE_VCF = file_vcf
    _WORKER_INDIVIDUAL_NAMES = individual_names


def processOneLineOfVCFFast(line):
    es = line.rstrip('\n').split('\t')
    if len(es) < 8:
        return processOneLineOfVCF(line, _WORKER_INDIVIDUAL_COL, _WORKER_CHROMOSOME_ONLY, _WORKER_FILTER_PASS, _WORKER_INFO_FIELD, _WORKER_INFO_FIELD_THRES, CHROMOSOMES, _WORKER_FILE_VCF, _WORKER_INDIVIDUAL_NAMES)

    chromosome, position, reference, alternatives, FILTER = es[0], es[1], es[3], es[4], es[6]

    if _WORKER_CHROMOSOME_ONLY and chromosome not in CHROMOSOMES:
        return ''
    if _WORKER_FILTER_PASS and FILTER != 'PASS':
        return ''
    if _WORKER_INFO_FIELD and _WORKER_INFO_FIELD_THRES:
        passes_threshold = info_field_passes_threshold(es[7], _WORKER_INFO_FIELD, _WORKER_INFO_FIELD_THRES)
        if passes_threshold is None:
            print('waring! INFO field', _WORKER_INFO_FIELD, 'not found in', _WORKER_FILE_VCF, line)
        elif not passes_threshold:
            return ''

    if ',' in alternatives:
        return processOneLineOfVCF(line, _WORKER_INDIVIDUAL_COL, _WORKER_CHROMOSOME_ONLY, _WORKER_FILTER_PASS, _WORKER_INFO_FIELD, _WORKER_INFO_FIELD_THRES, CHROMOSOMES, _WORKER_FILE_VCF, _WORKER_INDIVIDUAL_NAMES)
    alternative = alternatives
    if alternative == '*':
        return ''
    if len(reference) > 1 and len(alternative) > 1:
        return processOneLineOfVCF(line, _WORKER_INDIVIDUAL_COL, _WORKER_CHROMOSOME_ONLY, _WORKER_FILTER_PASS, _WORKER_INFO_FIELD, _WORKER_INFO_FIELD_THRES, CHROMOSOMES, _WORKER_FILE_VCF, _WORKER_INDIVIDUAL_NAMES)

    if not _WORKER_INDIVIDUAL_COL:
        return f'{chromosome}\t{position}\t{reference}\t{alternative}\n'

    genotypes = _decode_record_genotypes(es, _WORKER_INDIVIDUAL_COL, _WORKER_FILE_VCF, _WORKER_INDIVIDUAL_NAMES)
    if all(all(allele == 0 for allele in genotype) for genotype in genotypes):
        return ''
    bits = [
        str(1 if allele == 1 else 0)
        for genotype in genotypes
        for allele in genotype
    ]
    return f'{chromosome}\t{position}\t{reference}\t{alternative}\t' + '\t'.join(bits) + '\n'


def processVCFChunkFast(lines):
    return ''.join(processOneLineOfVCFFast(line) for line in lines)


def vcf_line_chunks(file_vcf, chunk_size=2000, threads=1):
    chunk = []
    for line in iter_vcf_records(file_vcf, threads=threads):
        chunk.append(line)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def convertVCF2MutationComplexStream(file_vcf, fout, individual_col, chromosome_only=True, filter_PASS=True, info_field=None, info_field_thres=None, threads=1, chunk_size=2000, individual_names=None):
    '''Convert VCF records to mutation TSV through an ordered streaming worker pool.'''
    if threads <= 1:
        init_vcf_worker(individual_col, chromosome_only, filter_PASS, info_field, info_field_thres, file_vcf, individual_names)
        for lines in vcf_line_chunks(file_vcf, chunk_size=chunk_size, threads=threads):
            fout.write(processVCFChunkFast(lines))
        return

    pool = Pool(threads, initializer=init_vcf_worker, initargs=(individual_col, chromosome_only, filter_PASS, info_field, info_field_thres, file_vcf, individual_names))
    chunks = vcf_line_chunks(file_vcf, chunk_size=chunk_size, threads=threads)
    try:
        results = pool.imap(processVCFChunkFast, chunks, chunksize=1)
        if tqdm:
            results = tqdm.tqdm(results, desc=f"vcf converting to tsv: {os.path.basename(file_vcf)}")
        for text in results:
            if text:
                fout.write(text)
    finally:
        pool.close()
        pool.join()


def process_manifest_row_python(task):
    '''Convert one manifest row/VCF into variant keys and two allele-presence bits.'''
    row_index, row, chromosome_only, filter_PASS, info_field, info_field_thres = task
    individual_col = [row['sample_col']]
    row_values = {}
    for line in iter_vcf_records(row['filepath'], threads=1):
        new_lines = processOneLineOfVCF(line, individual_col, chromosome_only, filter_PASS, info_field, info_field_thres, CHROMOSOMES, row['filepath'], [row['sample']])
        if not new_lines:
            continue
        for new_line in new_lines.rstrip('\n').split('\n'):
            fields = new_line.split('\t')
            key = tuple(fields[:4])
            values = fields[4:]
            if key not in row_values:
                row_values[key] = ['0', '0']
            for offset, value in enumerate(values[:2]):
                if value == '1':
                    row_values[key][offset] = '1'
    return row_index, row['filepath'], row_values


def merge_manifest_row_values(variant_to_values, row_index, row_values, n_columns):
    '''Merge one manifest worker result into the population allele matrix.'''
    column_offset = row_index * 2
    for key, values in row_values.items():
        if key not in variant_to_values:
            variant_to_values[key] = ['0'] * n_columns
        for offset, value in enumerate(values[:2]):
            if value == '1':
                variant_to_values[key][column_offset + offset] = '1'


def convert_manifest_with_python(rows, file_output, file_output_done, filter_PASS=True, chromosome_only=True, info_field=None, info_field_thres=None, threads=1):
    '''Convert manifest VCF rows directly into a population mutation TSV matrix.'''
    column_for_samples = []
    for row in rows:
        column_for_samples.append(row['name_use'] + '__1')
        column_for_samples.append(row['name_use'] + '__2')
    variant_to_values = {}
    n_columns = len(column_for_samples)
    tasks = [(row_index, row, chromosome_only, filter_PASS, info_field, info_field_thres) for row_index, row in enumerate(rows)]
    if threads > 1 and len(tasks) > 1:
        pool = Pool(min(threads, len(tasks)))
        try:
            results = pool.imap_unordered(process_manifest_row_python, tasks, chunksize=1)
            if tqdm:
                results = tqdm.tqdm(results, total=len(tasks), desc='manifest vcf converting to tsv')
            for row_index, file_vcf, row_values in results:
                merge_manifest_row_values(variant_to_values, row_index, row_values, n_columns)
                print('finish converting vcf file:', file_vcf)
        finally:
            pool.close()
            pool.join()
    else:
        to_iter = tasks
        if tqdm:
            to_iter = tqdm.tqdm(to_iter, total=len(tasks), desc='manifest vcf converting to tsv')
        for task in to_iter:
            row_index, file_vcf, row_values = process_manifest_row_python(task)
            merge_manifest_row_values(variant_to_values, row_index, row_values, n_columns)
            print('finish converting vcf file:', file_vcf)

    with open(file_output, 'w') as fout:
        fout.write('\t'.join(['chr', 'pos', 'ref', 'alt'] + column_for_samples) + '\n')
        for key in sorted(variant_to_values, key=lambda x: (x[0], int(x[1]) if str(x[1]).isdigit() else x[1], x[2], x[3])):
            fout.write('\t'.join(list(key) + variant_to_values[key]) + '\n')
    open(file_output_done, 'w').write('\n'.join(column_for_samples))
    return column_for_samples


def convertVCFManifest2MutationComplex(file_manifest, outprefix=None, filter_PASS=True, chromosome_only=True, info_field=None, info_field_thres=None, threads=1, force=False):
    '''Convert a filepath manifest into the same population TSV produced from multi-sample VCF input.'''
    rows = read_vcf_manifest(file_manifest)
    if outprefix is None:
        file_output = strip_vcf_suffix(file_manifest) + '.tsv'
        outprefix = os.path.splitext(file_output)[0]
        print(f"Output prefix not specified, writing to '{file_output}'")
    else:
        file_output = outprefix + '.tsv'
    file_output_done = file_output + '.done'
    cache = StageCache(
        file_output, [file_manifest] + [row['filepath'] for row in rows],
        {'mode': 'manifest', 'rows': rows, 'filter_pass': filter_PASS,
         'chromosome_only': chromosome_only, 'info_field': info_field,
         'info_field_thres': info_field_thres, 'threads': threads},
        force=force,
    )
    if cache.prepare():
        print(f"Output file '{file_output}' already exists. Skipping.")
        return open(file_output_done, 'r').read().strip().split()

    if info_field:
        try:
            info_field_thres = float(info_field_thres)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"--info_field_thres ({info_field_thres!r}) must be a number"
            ) from exc

    print('use Python parser for manifest VCF files.')
    columns = convert_manifest_with_python(rows, file_output, file_output_done, filter_PASS, chromosome_only, info_field, info_field_thres, threads)
    cache.complete()
    return columns



def tsv2memmap(tsv_file, individuals = None, memmap_file=None, batch_size=100, force=False):
    '''tsv file is a tsv file with columns: chr, pos, ref, sample1__1, sample1__2, sample2__1, sample2__2, ..., sampleN__1, sampleN__2 columns
    individuals is a list of individuals to be used. the column values must be 0 or 1
    if memmap_file is None, memmap_file = tsv_file + '.memmap'
    '''
    if memmap_file is None:
        memmap_file = tsv_file + '.memmap'
    
    memmap_file_done = memmap_file + '.done'

    with open(tsv_file, 'rb') as fo:
        header = fo.readline().decode().rstrip('\n').split('\t')

    if individuals is None:
        individuals = [i for i in header if i not in ['chr', 'pos', 'ref', 'alt','pos_end']]
    elif isinstance(individuals, str) and ',' in individuals:
        individuals = individuals.split(',')
    elif isinstance(individuals, str):
        individuals = [individuals]
    else:
        individuals = individuals

    column_to_index = {column: i for i, column in enumerate(header)}
    missing = [individual for individual in individuals if individual not in column_to_index]
    if missing:
        raise ValueError('individual columns not found in tsv file: ' + ','.join(missing))
    if not individuals:
        raise ValueError('no individual columns selected for memmap conversion')
    cache = StageCache(
        memmap_file, [tsv_file],
        {'individuals': list(individuals), 'dtype': 'int8'}, force=force,
    )
    if cache.prepare():
        print(f"memmap file '{memmap_file}' already exists, skipping")
        return memmap_file
    individual_indices = [column_to_index[individual] for individual in individuals]
    is_contiguous = individual_indices == list(range(individual_indices[0], individual_indices[0] + len(individual_indices)))

    with open(tsv_file, 'rb') as fo:
        total_rows = sum(block.count(b'\n') for block in iter(lambda: fo.read(8 * 1024 * 1024), b'')) - 1

    shape = (total_rows, len(individuals))
    mmap = np.memmap(memmap_file, dtype='int8', mode='w+', shape=shape)

    print('storing tsv file to memmap...')

    def field_start(line, column_index):
        if column_index == 0:
            return 0
        pos = -1
        for _ in range(column_index):
            pos = line.find(b'\t', pos + 1)
            if pos == -1:
                return -1
        return pos + 1

    with open(tsv_file, 'rb') as fo:
        fo.readline()
        to_iter = enumerate(fo)
        if tqdm:
            to_iter = tqdm.tqdm(to_iter, total=total_rows)
        for row_offset, line in to_iter:
            line = line.rstrip(b'\n\r')
            if is_contiguous:
                start = field_start(line, individual_indices[0])
                if start == -1:
                    raise ValueError('malformed TSV line with too few columns')
                if individual_indices[-1] == len(header) - 1:
                    segment = line[start:]
                else:
                    end = field_start(line, individual_indices[-1] + 1)
                    if end == -1:
                        raise ValueError('malformed TSV line with too few columns')
                    segment = line[start:end - 1]
                values = np.frombuffer(segment, dtype=np.uint8)[::2] - ord('0')
                if values.shape[0] != len(individuals):
                    fields = line.split(b'\t')
                    values = np.fromiter((int(fields[i]) for i in individual_indices), dtype=np.int8, count=len(individual_indices))
                mmap[row_offset] = values
            else:
                fields = line.split(b'\t')
                mmap[row_offset] = np.fromiter((int(fields[i]) for i in individual_indices), dtype=np.int8, count=len(individual_indices))

    mmap.flush()
    del mmap
    print(f"TSV file '{tsv_file}' has been successfully converted to memory-mapped file '{memmap_file}'")
    open(memmap_file_done, 'w').close()
    cache.complete()
    return memmap_file

def get_header_sample_col(file_vcf,individual_input="ALL_SAMPLES"):
    '''
    return the header line based on files_vcf and individual_input
    file_vcf should be a vcf file
    return header line and individual_col
    '''
    vcf_file_path = file_vcf
    try:
        fo = openFile(vcf_file_path)
    except FileNotFoundError:
        print(f"Warning: Could not open file {vcf_file_path}. Skipping.")
        return '', [], None
    except Exception as e:
            print(f"Warning: Error opening file {vcf_file_path}: {e}. Skipping.")
            return '', [], None
    
    # Read header lines to find samples
    header_line = None
    for line in fo:
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            header_line = line
            break
        else: # Unexpected line before #CHROM header
            print(f"Warning: Unexpected line format before #CHROM header in {vcf_file_path}. Attempting to continue.")
            header_line = line # Try processing it as header anyway? Risky.
            break

    if header_line is None:
        print(f"Warning: #CHROM header line not found in {vcf_file_path}. Skipping.")
        fo.close()
        return '', [], None
            
    columns = header_line.strip().split('\t')
    if len(columns) < 8:
        print(f"Warning: #CHROM header line in {vcf_file_path} has fewer than 8 columns. Skipping.")
        fo.close()
        return '', [], None
    try:
        individual, individual_col = _resolve_sample_selection(columns, individual_input, vcf_file_path)
    except Exception:
        fo.close()
        raise
    
    
    # Define output column names (only needs to be done once if consistent across files)
    current_column_for_samples = []
    for indi_name in individual: # Use determined 'individual' list
        current_column_for_samples.append(indi_name + '__1')
        current_column_for_samples.append(indi_name + '__2')

    output_column_names = ['chr', 'pos', 'ref', 'alt'] + current_column_for_samples
    header_string = '\t'.join(output_column_names) + '\n'
    
    return header_string, individual_col, fo

def convertVCF2MutationComplex(file_vcf, outprefix = None, individual_input="ALL_SAMPLES", filter_PASS = True, chromosome_only = True, info_field = None, info_field_thres=None, threads = 1, force=False):
    '''convert vcf file to tsv file. with columns: chr, pos, ref, sample1__1, sample1__2, sample2__1, sample2__2, ..., sampleN__1, sampleN__2 columns
    If individual is None, use the first sample in the vcf file. 
    If individual == 'ALL_SAMPLES', use all samples in the vcf file.
    If individual == 'ALL_VARIANTS', ignore sample columns
    info_field is the INFO field key used to filter. In ProHap, the default setting is "AF"
    info_field_thres is the threshold of the INFO field key used to filter. In ProHap, the default setting is 0.01
    
    '''
    if is_manifest_file(file_vcf):
        return convertVCFManifest2MutationComplex(file_vcf, outprefix, filter_PASS, chromosome_only, info_field, info_field_thres, threads, force=force)

    if ',' in file_vcf:
        files_vcf = file_vcf.split(',')
    elif os.path.exists(file_vcf):
        files_vcf = [file_vcf]
    else:
        files_vcf = glob.glob(file_vcf)
        if len(files_vcf) == 0:
            print(file_vcf, 'not found.')
            return []
    
    if outprefix is None:
        # Derive output name from first VCF file if multiple, handle potential .gz
        base_name = os.path.basename(files_vcf[0])
        if base_name.endswith(".vcf.gz"):
             file_output = base_name[:-7] + '.tsv'
        elif base_name.endswith(".vcf"):
             file_output = base_name[:-4] + '.tsv'
        else:
             file_output = base_name + '.tsv'
        print(f"Output prefix not specified, writing to '{file_output}'")
        # file_output = file_vcf + '.tsv'
    else:
        file_output = outprefix + '.tsv'
    
    # check if file_output is already finished
    file_output_done = file_output + '.done'
    cache = StageCache(
        file_output, files_vcf,
        {'mode': 'vcf', 'selection': individual_input,
         'filter_pass': filter_PASS, 'chromosome_only': chromosome_only,
         'info_field': info_field, 'info_field_thres': info_field_thres,
         'threads': threads},
        force=force,
    )
    if cache.prepare():
        print(f"Output file '{file_output}' already exists. Skipping.")
        return open(file_output_done, 'r').read().strip().split()
    
    if info_field:
        try:
            info_field_thres = float(info_field_thres)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"--info_field_thres ({info_field_thres!r}) must be a number"
            ) from exc

    dc_header_string_and_other_info = {
        file_vcf: get_header_sample_col(file_vcf, individual_input=individual_input)
        for file_vcf in files_vcf
    }
    dc_header_string_and_other_info = {k:v for k,v in dc_header_string_and_other_info.items() if v[0]}
    if len(dc_header_string_and_other_info) == 0:
        return []
    if len(set([i[0] for i in dc_header_string_and_other_info.values()])) != 1:
        print(f"Error: Different header strings found in {files_vcf}.")
        for _, _, fo in dc_header_string_and_other_info.values():
            fo.close()
        return []
    header_string = list(dc_header_string_and_other_info.values())[0][0]
    column_for_samples = header_string.strip().split()[4:]
    individual_names = [column[:-3] for column in column_for_samples[::2]]
    
    total_vcf_size = sum([os.path.getsize(file_vcf) for file_vcf in files_vcf])
    use_stream_pool = threads > 1 and total_vcf_size >= 10*1024*1024

    with open(file_output, 'w') as fout:
        fout.write(header_string)
        if use_stream_pool:
            print('use multiple threads to convert the vcf files by streaming chunks')
            for file_vcf in dc_header_string_and_other_info:
                _, individual_col, fo = dc_header_string_and_other_info[file_vcf]
                fo.close()
                convertVCF2MutationComplexStream(file_vcf, fout, individual_col, chromosome_only, filter_PASS, info_field, info_field_thres, threads=threads, individual_names=individual_names)
                print('finish converting vcf file:', file_vcf)
        else:
            for file_vcf in dc_header_string_and_other_info:
                _, individual_col, fo = dc_header_string_and_other_info[file_vcf]
                if tqdm:
                    to_iter = tqdm.tqdm(fo)
                else:
                    to_iter = fo
                for line in to_iter:
                    new_line = processOneLineOfVCF(line, individual_col, chromosome_only, filter_PASS, info_field, info_field_thres, CHROMOSOMES, file_vcf, individual_names)
                    if new_line:
                        fout.write(new_line)
                fo.close()
                print('finish converting vcf file:', file_vcf)

    open(file_output_done, 'w').write('\n'.join(column_for_samples))
    cache.complete()
    return column_for_samples

description = '''convert extract mutation information from vcf file
If --sample is not set or only one sample were used, the same as PrecisionProDB version 1.0. two files will be generated.
Otherwise, PrecisionProDB version 2.0 mode. Output a single file with extra columns for each sample. If only one sample is used and version 2.0 mode is desired, use "--sample sample,"
'''

def build_parser():
    import argparse
    parser = argparse.ArgumentParser(description=description)
    add_argument_set(parser, 'vcf_conversion')
    return parser


def run_from_args(f):
    chromosome_only = not f.all_chromosomes
    filter_PASS = not f.no_filter
    sample = f.sample

    if is_manifest_file(f.file_vcf):
        print('convert manifest VCF list to mutation information in version 2.0 mode')
        convertVCF2MutationComplex(file_vcf = f.file_vcf, outprefix = f.outprefix, individual_input='ALL_SAMPLES', filter_PASS = filter_PASS, chromosome_only = chromosome_only, info_field = f.info_field, info_field_thres = f.info_field_thres, threads = f.threads, force=f.force)
    elif sample is None:
        print('convert vcf to mutation information in version 1.0 mode')
        getMutationsFromVCF(file_vcf = f.file_vcf, outprefix = f.outprefix, individual=f.sample, filter_PASS = filter_PASS, chromosome_only = chromosome_only)
    elif ',' in  sample or sample == 'ALL_SAMPLES' or sample == 'ALL_VARIANTS':
        print('convert vcf to mutation information in version 2.0 mode')
        convertVCF2MutationComplex(file_vcf = f.file_vcf, outprefix = f.outprefix, individual_input=f.sample, filter_PASS = filter_PASS, chromosome_only = chromosome_only, info_field = f.info_field, info_field_thres = f.info_field_thres, threads = f.threads, force=f.force)
    else:
        print('convert vcf to mutation information in version 1.0 mode')
        getMutationsFromVCF(file_vcf = f.file_vcf, outprefix = f.outprefix, individual=f.sample, filter_PASS = filter_PASS, chromosome_only = chromosome_only)


def main(argv=None):
    run_from_args(build_parser().parse_args(argv))


if __name__ == '__main__':
    main()
