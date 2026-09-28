import csv

import pytest

from precisionprodb.vcf2mutation import (
    convertVCF2MutationComplex,
    convertVCFManifest2MutationComplex,
    getMutationsFromVCF,
    init_vcf_worker,
    processOneLineOfVCF,
    processOneLineOfVCFFast,
)


def write_vcf(path, samples, records, format_fields='GT'):
    columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    if samples is not None:
        columns.extend(['FORMAT'] + samples)
    with open(path, 'w') as handle:
        handle.write('##fileformat=VCFv4.2\n')
        handle.write('\t'.join(columns) + '\n')
        for record in records:
            handle.write('\t'.join(record) + '\n')
    return path


def test_selected_sample_and_all_variants_are_honored(tmp_path):
    vcf = write_vcf(
        tmp_path / 'two_samples.vcf',
        ['sampleA', 'sampleB'],
        [
            ['1', '10', '.', 'A', 'G', '.', 'PASS', '.', 'GT', '0/0', '0/1'],
            ['1', '20', '.', 'C', 'T', '.', 'PASS', '.', 'GT', '0/0', '0/0'],
        ],
    )

    output_a = convertVCF2MutationComplex(str(vcf), str(tmp_path / 'sample_a'), individual_input='sampleA')
    rows_a = list(csv.reader(open(tmp_path / 'sample_a.tsv'), delimiter='\t'))
    assert output_a == ['sampleA__1', 'sampleA__2']
    assert rows_a == [['chr', 'pos', 'ref', 'alt', 'sampleA__1', 'sampleA__2']]

    convertVCF2MutationComplex(str(vcf), str(tmp_path / 'sample_b'), individual_input='sampleB')
    rows_b = list(csv.reader(open(tmp_path / 'sample_b.tsv'), delimiter='\t'))
    assert rows_b[1] == ['1', '10', 'A', 'G', '0', '1']

    columns = convertVCF2MutationComplex(
        str(vcf), str(tmp_path / 'all_variants'), individual_input='ALL_VARIANTS'
    )
    rows_all = list(csv.reader(open(tmp_path / 'all_variants.tsv'), delimiter='\t'))
    assert columns == []
    assert rows_all == [
        ['chr', 'pos', 'ref', 'alt'],
        ['1', '10', 'A', 'G'],
        ['1', '20', 'C', 'T'],
    ]


def test_explicit_missing_sample_fails_instead_of_using_all_variants(tmp_path):
    vcf = write_vcf(
        tmp_path / 'one_sample.vcf',
        ['present'],
        [['1', '10', '.', 'A', 'G', '.', 'PASS', '.', 'GT', '0/1']],
    )
    with pytest.raises(ValueError, match='not found.*absent'):
        convertVCF2MutationComplex(str(vcf), str(tmp_path / 'missing'), individual_input='absent')
    assert not (tmp_path / 'missing.tsv').exists()


@pytest.mark.parametrize(
    ('genotype', 'expected'),
    [
        ('1/.', '1\t0'),
        ('./1', '0\t1'),
        ('1', '1\t0'),
    ],
)
def test_missing_and_haploid_alleles_fill_reference_slot(genotype, expected):
    line = f'1\t10\t.\tA\tG\t.\tPASS\t.\tDP:GT\t7:{genotype}\n'
    assert processOneLineOfVCF(line, [9], file_vcf='fixture.vcf', individual_names=['sampleA']) == (
        f'1\t10\tA\tG\t{expected}\n'
    )


@pytest.mark.parametrize('genotype', ['.', './.'])
def test_all_missing_calls_are_reference_and_do_not_emit_variant(genotype):
    line = f'1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t{genotype}\n'
    assert processOneLineOfVCF(line, [9], file_vcf='fixture.vcf', individual_names=['sampleA']) == ''


def test_more_than_two_alleles_raise_with_sample_and_position_context():
    line = '1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1/1\n'
    with pytest.raises(ValueError, match='More than two.*sample sampleA.*fixture.vcf:1:10'):
        processOneLineOfVCF(line, [9], file_vcf='fixture.vcf', individual_names=['sampleA'])


def test_standard_and_fast_parsers_share_genotype_rules():
    line = '1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t1/.\n'
    expected = processOneLineOfVCF(line, [9], file_vcf='fixture.vcf', individual_names=['sampleA'])
    init_vcf_worker([9], True, True, None, None, 'fixture.vcf', ['sampleA'])
    assert processOneLineOfVCFFast(line) == expected


def test_manifest_parser_uses_same_haploid_padding(tmp_path):
    vcf = write_vcf(
        tmp_path / 'haploid.vcf',
        ['source_sample'],
        [['1', '10', '.', 'A', 'G', '.', 'PASS', '.', 'GT', '1']],
    )
    manifest = tmp_path / 'manifest.tsv'
    manifest.write_text(f'filepath\tsample\tname_use\n{vcf}\tsource_sample\tcohortA\n')
    convertVCFManifest2MutationComplex(str(manifest), str(tmp_path / 'manifest_out'))
    rows = list(csv.reader(open(tmp_path / 'manifest_out.tsv'), delimiter='\t'))
    assert rows[0] == ['chr', 'pos', 'ref', 'alt', 'cohortA__1', 'cohortA__2']
    assert rows[1] == ['1', '10', 'A', 'G', '1', '0']


def test_legacy_parser_respects_sample_and_genotype_policy(tmp_path):
    vcf = write_vcf(
        tmp_path / 'legacy.vcf',
        ['sampleA', 'sampleB'],
        [['1', '10', '.', 'A', 'G', '.', 'PASS', '.', 'GT', '0/0', '1/.']],
    )
    df1, df2 = getMutationsFromVCF(str(vcf), individual='sampleB')
    assert df1[['pos', 'alt']].values.tolist() == [[10, 'G']]
    assert df2.empty


def test_header_only_vcf_without_sample_columns_can_convert_all_variants(tmp_path):
    vcf = write_vcf(
        tmp_path / 'sites_only.vcf',
        None,
        [['1', '10', '.', 'A', 'G', '.', 'PASS', '.']],
    )
    convertVCF2MutationComplex(str(vcf), str(tmp_path / 'sites_only_out'), individual_input='ALL_VARIANTS')
    rows = list(csv.reader(open(tmp_path / 'sites_only_out.tsv'), delimiter='\t'))
    assert rows == [['chr', 'pos', 'ref', 'alt'], ['1', '10', 'A', 'G']]
