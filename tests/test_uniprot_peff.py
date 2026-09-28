"""UniProt PEFF carries annotations from every matching reference ID."""

from precisionprodb.generatePEFFoutput import generateUniprotPEFFout


def test_uniprot_peff_combines_tied_reference_variants(tmp_path):
    uniprot = tmp_path / 'uniprot.fa'
    uniprot.write_text('>U\nMPEPTIDEK\n')
    references = tmp_path / 'refs.peff.fa'
    references.write_text(
        '# PEFF 1.0\n'
        '>PrecisionProDB:P1 \\VariantSimple=(2|V|1-4-A-G) \\Length=9\nMPEPTIDEK\n'
        '>PrecisionProDB:P2 \\VariantSimple=(5|A|1-8-C-T) \\Length=9\nMPEPTIDEK\n'
    )
    relationships = tmp_path / 'changed.tsv'
    relationships.write_text('uniprot_id\tref_id\nU\tP2\nU\tP1\n')
    output = tmp_path / 'out.peff.fa'
    generateUniprotPEFFout(
        str(references), str(uniprot), str(relationships), str(output)
    )
    result = output.read_text()
    assert '>PrecisionProDB:U \\VariantSimple=(2|V|1-4-A-G)(5|A|1-8-C-T)' in result
    assert '\\Length=9\nMPEPTIDEK\n' in result
