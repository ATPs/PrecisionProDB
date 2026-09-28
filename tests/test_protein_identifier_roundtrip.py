"""Protein identifiers must survive annotation TSV round trips verbatim."""

import csv
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO

from precisionprodb.PrecisionProDB import main
from precisionprodb import PrecisionProDB_vcf as vcf_pipeline
from precisionprodb import perChrom
from precisionprodb.PrecisionProDB_core import PerGeno
from precisionprodb.generatePEFFoutput import generatePEFFoutput


@pytest.fixture(params=['123', '00123', 'NA'])
def protein_id(request):
    return request.param


def save_changed(prefix, protein_id):
    rows = pd.DataFrame(
        [{
            'protein_id_fasta': protein_id,
            'protein_description': protein_id + ' description',
            'AA_seq': 'MAK',
            'new_AA': 'MVK',
            'variant_AA': 'A2V(1-5-C-T)',
            'insertion_AA': None,
        }],
        index=pd.Index([protein_id], name='protein_id'),
    )
    perChrom.save_mutation_and_proteins(rows, str(prefix))


def read_annotations(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def test_chromosome_merge_preserves_literal_ids(tmp_path, monkeypatch, protein_id):
    tempfolder = tmp_path / 'out_temp'
    tempfolder.mkdir()
    save_changed(tempfolder / '1', protein_id)
    (tempfolder / '1.proteins.fa').write_text(
        f'>{protein_id}\t{protein_id} description\nMAK\n'
    )
    (tempfolder / 'not_in_genome.proteins.fa').write_text('')
    pipeline = PerGeno.__new__(PerGeno)
    pipeline.outprefix = str(tmp_path / 'out')
    pipeline.tempfolder = str(tempfolder)
    pipeline.chromosomes = ['1']
    pipeline.keep_all = True
    monkeypatch.setattr(pipeline, 'runSinglePerChrom', lambda chromosome: chromosome)

    pipeline.runPerChom()

    rows = read_annotations(tmp_path / 'out.pergeno.aa_mutations.csv')
    assert rows[0]['protein_id_fasta'] == protein_id
    assert rows[0]['protein_id'] == protein_id
    assert rows[0]['insertion_AA'] == ''
    records = list(SeqIO.parse(tmp_path / 'out.pergeno.protein_changed.fa', 'fasta'))
    assert [(record.id, str(record.seq)) for record in records] == [(protein_id, 'MVK')]


def test_legacy_vcf_merge_preserves_literal_ids(tmp_path, monkeypatch, protein_id):
    class PreparedHaplotype:
        def __init__(self, **kwargs):
            self.prefix = kwargs['outprefix']

        def splitInputByChromosomes(self):
            pass

        def runPerChom(self):
            prefix = self.prefix + '.pergeno'
            save_changed(prefix, protein_id)
            Path(prefix + '.mutated_protein.fa').rename(
                prefix + '.protein_changed.fa'
            )

    monkeypatch.setattr(vcf_pipeline, 'PerGeno', PreparedHaplotype)
    monkeypatch.setattr(vcf_pipeline, 'getMutationsFromVCF', lambda **kwargs: None)
    reference = tmp_path / 'reference.fa'
    reference.write_text(f'>{protein_id} description\nMAK\n')

    vcf_pipeline.runPerGenoVCF(
        '', '', '', str(reference), threads=1, outprefix=str(tmp_path / 'out'),
    )

    rows = read_annotations(tmp_path / 'out.pergeno.aa_mutations.csv')
    expected_id = protein_id + '__12'
    assert [row['protein_id_fasta'] for row in rows] == [expected_id]
    assert rows[0]['protein_id'] == protein_id
    for filename in ('out.pergeno.protein_changed.fa', 'out.pergeno.protein_all.fa'):
        records = list(SeqIO.parse(tmp_path / filename, 'fasta'))
        assert [(record.id, str(record.seq)) for record in records] == [(expected_id, 'MVK')]


def test_peff_preserves_literal_ids_and_variant_annotations(tmp_path, protein_id):
    prefix = tmp_path / 'out'
    save_changed(prefix, protein_id)
    reference = tmp_path / 'reference.fa'
    reference.write_text(f'>{protein_id} description\nMAK\n')
    output = tmp_path / 'out.peff.fa'

    generatePEFFoutput(str(reference), str(prefix) + '.aa_mutations.csv', str(output))

    assert (
        f'>PrecisionProDB:{protein_id} \\VariantSimple=(2|V|1-5-C-T)'
        ' \\Length=3\nMAK\n'
    ) in output.read_text()


def test_numeric_and_literal_ids_through_legacy_vcf_cli(tmp_path, protein_id):
    genome = tmp_path / 'genome.fa'
    genome.write_text('>1\nATGGCTAAATAA\n')
    protein = tmp_path / 'protein.fa'
    protein.write_text(f'>{protein_id}\nMAK\n')
    gtf = tmp_path / 'genes.gtf'
    attributes = f'gene_id "gene1"; transcript_id "{protein_id}";'
    gtf.write_text(
        f'1\tfixture\texon\t1\t12\t.\t+\t.\t{attributes}\n'
        f'1\tfixture\tCDS\t1\t9\t.\t+\t0\t{attributes}\n'
    )
    vcf = tmp_path / 'variants.vcf'
    vcf.write_text(
        '##fileformat=VCFv4.2\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n'
        '1\t5\t.\tC\tT\t.\tPASS\t.\tGT\t1|1\n'
    )

    main([
        '-g', str(genome), '-p', str(protein), '-f', str(gtf),
        '-m', str(vcf), '-o', str(tmp_path / 'out'), '-a', 'gtf',
        '-k', 'transcript_id', '--sqlite', 'NONE', '--PEFF', '-t', '1',
    ])

    expected_id = protein_id + '__12'
    rows = read_annotations(tmp_path / 'out.pergeno.aa_mutations.csv')
    assert [row['protein_id_fasta'] for row in rows] == [expected_id]
    records = list(SeqIO.parse(tmp_path / 'out.pergeno.protein_changed.fa', 'fasta'))
    assert [(record.id, str(record.seq)) for record in records] == [(expected_id, 'MVK')]
    assert (
        f'>PrecisionProDB:{protein_id} \\VariantSimple=(2|V|1-5-C-T)'
    ) in (tmp_path / 'out.pergeno.protein_PEFF.fa').read_text()
