import argparse
import csv
import os
import shlex
import subprocess
import sys
import time

from Bio import SeqIO


def run_command(command, cwd=None):
    """
    Run a system command using subprocess and print the output.
    """
    display_command = ' '.join(shlex.quote(part) for part in command)
    print(f"Running command: {display_command}")
    start_time = time.time()
    result = subprocess.run(command, cwd=cwd, check=False)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Command finished in {elapsed_time:.2f} seconds\n\n")
    if result.returncode != 0:
        print(f"Error running command (exit {result.returncode}): {display_command}\n\n")
        return display_command, 'failed'
    else:
        return display_command, 'completed'


def validate_case_outputs(output_prefix, key_input, expect_uniprot=False):
    """Check core artifacts and their identifier relationships."""
    mutation_file = output_prefix + '.pergeno.aa_mutations.csv'
    changed_file = output_prefix + '.pergeno.protein_changed.fa'
    all_file = output_prefix + '.pergeno.protein_all.fa'
    peff_file = output_prefix + '.pergeno.protein_PEFF.fa'
    required = [mutation_file, changed_file, all_file, peff_file]
    if expect_uniprot:
        required.extend([
            output_prefix + '.uniprot_changed.tsv',
            output_prefix + '.uniprot_changed.fa',
            output_prefix + '.uniprot_all.fa',
            output_prefix + '.uniprot_PEFF.fa',
        ])
    missing = [path for path in required if not os.path.isfile(path)]
    if missing:
        raise RuntimeError('missing output files: ' + ', '.join(missing))

    changed_records = list(SeqIO.parse(changed_file, 'fasta'))
    all_records = list(SeqIO.parse(all_file, 'fasta'))
    changed_entries = {(record.id, str(record.seq)) for record in changed_records}
    all_entries = {(record.id, str(record.seq)) for record in all_records}
    if not changed_entries.issubset(all_entries):
        raise RuntimeError('changed protein FASTA contains records absent from all-protein FASTA')

    with open(mutation_file, newline='') as handle:
        rows = csv.DictReader(handle, delimiter='\t')
        if not rows.fieldnames or 'protein_id_fasta' not in rows.fieldnames:
            raise RuntimeError('mutation table is missing the protein_id_fasta column')
        mutated_ids = {row['protein_id_fasta'] for row in rows if row.get('protein_id_fasta')}
    changed_ids = {record.id for record in changed_records}
    absent = sorted(mutated_ids - changed_ids)
    if absent:
        raise RuntimeError('mutation table proteins absent from changed FASTA: ' + ', '.join(absent[:10]))

    # Parse PEFF and UniProt FASTAs as well; empty sequence files are valid.
    for path in [peff_file] + required[4:]:
        if path.endswith('.fa'):
            list(SeqIO.parse(path, 'fasta'))


def get_cmd_to_run(key_input, key_variant, sqlite_key, output_test, dc_variant, dc_inputs, path_of_precisionprodb):
    '''
    '''
    folder_work = os.path.join(output_test, key_input, key_variant, sqlite_key)
    if not os.path.exists(folder_work):
        os.makedirs(folder_work)
    file_mutation = dc_variant[key_variant]
    file_genome = dc_inputs[key_input]['genome']
    file_protein = dc_inputs[key_input]['protein']
    file_gtf = dc_inputs[key_input]['gtf']
    datatype = dc_inputs[key_input]['datatype']
    output = os.path.join(folder_work, f'{key_input}.{key_variant}.{sqlite_key}')
    file_sqlite = output + '.sqlite'
    script_folder = os.path.join(path_of_precisionprodb, 'src', 'precisionprodb')
    
    print(f'running test with key_input: {key_input}, key_variant: {key_variant}, sqlite_key: {sqlite_key}')
    precisionprodb_script = os.path.join(script_folder, 'PrecisionProDB.py')
    python = sys.executable
    common = [
        python, precisionprodb_script,
        '-m', file_mutation,
        '-g', file_genome,
        '-p', file_protein,
        '-f', file_gtf,
        '-o', output,
        '-a', datatype,
        '--PEFF',
    ]
    if key_input == 'UniProt' and key_variant != 'str':
        common.extend(['-U', dc_inputs[key_input]['UniProt'], '-t', '4', '-D', 'Uniprot'])
    else:
        common.extend(['-t', '4'])

    if sqlite_key == 'no_sqlite':
        print("Running test: without use sqlite file")
        # SQLite is enabled by default; NONE selects the legacy path.
        return [[*common, '--sqlite', 'NONE']]
    elif sqlite_key == 'sqlite_one_step':
        print("Running test: use SQLite file as intermediate file")
        return [[*common, '-S', file_sqlite]]
    elif sqlite_key == 'sqlite_two_step':
        print("Running test: Generate SQLite file in advance and use SQLite")
        build_script = os.path.join(script_folder, 'buildSqlite.py')
        build_command = [
            python, build_script,
            '-S', file_sqlite,
            '-o', output + '.build',
            '-g', file_genome,
            '-p', file_protein,
            '-f', file_gtf,
            '-a', datatype,
        ]
        execute_command = [*common, '-S', file_sqlite]
        return [build_command, execute_command]

    raise ValueError(f'Unknown SQLite test mode: {sqlite_key}')

def main_test(dc_variant, dc_inputs, dc_sqlite, output_test, path_of_precisionprodb):

    ls_results = []
    failed = False
    for key_input in dc_inputs:
        for key_variant in dc_variant:
            for sqlite_key in dc_sqlite:
                if sqlite_key == 'no_sqlite' and key_variant == 'str':
                    continue
                commands = get_cmd_to_run(
                    key_input, key_variant, sqlite_key, output_test, dc_variant,
                    dc_inputs, path_of_precisionprodb,
                )
                folder_work = os.path.join(output_test, key_input, key_variant, sqlite_key)
                case_ok = True
                for command in commands:
                    result = run_command(command, cwd=folder_work)
                    ls_results.append(result)
                    if result[1] != 'completed':
                        case_ok = False
                        break
                if case_ok:
                    output_prefix = os.path.join(
                        folder_work, f'{key_input}.{key_variant}.{sqlite_key}'
                    )
                    try:
                        validate_case_outputs(
                            output_prefix,
                            key_input,
                            expect_uniprot=(key_input == 'UniProt' and key_variant != 'str'),
                        )
                        ls_results.append((output_prefix, 'outputs_validated'))
                    except Exception as error:
                        failed = True
                        print(f'Output validation failed for {output_prefix}: {error}')
                        ls_results.append((output_prefix, f'validation_failed: {error}'))
                else:
                    failed = True

    # save ls_results to file
    file_job_summary = os.path.join(output_test, 'test_running_summary.txt')
    with open(file_job_summary, 'w') as f:
        for i, result in enumerate(ls_results):
            f.write(f'{result[0]}\n{result[1]}\n\n')

    return 1 if failed else 0



def main():
    description = """
    Run PrecisionProDB tests.
    Usage: test.py -s PATH_OF_PRECISIONPRODB [-o OUTPUT_TEST]
    Make sure Python is in your PATH

    Options:
        -s PATH_OF_PRECISIONPRODB  Path to PrecisionProDB. Will use PATH_OF_PRECISIONPRODB/src/precisionprodb/ to find the scripts
                                and PATH_OF_PRECISIONPRODB/examples to find the test input files. if not set, will be determined based on the PrecisionProDB_test.py file. 
        -o OUTPUT_TEST             Path to store the output test results. If not set, will use PATH_OF_PRECISIONPRODB/test_output
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-s", "--src", required=False, help="Path to PrecisionProDB", default="")
    parser.add_argument("-o", "--output", help="Path to store output test results")

    TEST = False
    # TEST = True
    if TEST:
        args = parser.parse_args('-s /data/p/xiaolong/PrecisionProDB -o /XCLabServer002_fastIO/test_output/'.split())
    else:
        args = parser.parse_args()

    path_of_precisionprodb = args.src
    path_of_precisionprodb = (path_of_precisionprodb if path_of_precisionprodb else os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    path_of_precisionprodb = os.path.abspath(path_of_precisionprodb)
    output_test = (
        args.output if args.output else os.path.join(path_of_precisionprodb, "test_output")
    )
    if not os.path.isabs(output_test):
        output_test = os.path.abspath(output_test)

    print(f"PATH_OF_PRECISIONPRODB: {path_of_precisionprodb}")
    print(f"OUTPUT_TEST: {output_test}")

    # Define input files
    dc_variant = {}
    dc_variant["tsv"] = os.path.join(path_of_precisionprodb, "examples", "gnomAD.variant.txt.gz")
    dc_variant['vcf'] = os.path.join(path_of_precisionprodb, "examples", "celline.vcf.gz")
    dc_variant['str'] = "chr1-942451-T-C,1-6253878-C-T,1-2194700-C-G,1-1719406-G-A"

    dc_inputs = {}
    dc_inputs['Ensembl'] = {
        'gtf': os.path.join(path_of_precisionprodb, "examples", "Ensembl", "Ensembl.gtf.gz"),
        'genome': os.path.join(path_of_precisionprodb, "examples", "Ensembl", "Ensembl.genome.fa.gz"),
        'protein': os.path.join(path_of_precisionprodb, "examples", "Ensembl", "Ensembl.protein.fa.gz"),
        'datatype':'Ensembl_GTF'
    }
    dc_inputs['GENCODE'] = {
        'gtf': os.path.join(path_of_precisionprodb, "examples", "GENCODE", "GENCODE.gtf.gz"),
        'genome': os.path.join(path_of_precisionprodb, "examples", "GENCODE", "GENCODE.genome.fa.gz"),
        'protein': os.path.join(path_of_precisionprodb, "examples", "GENCODE", "GENCODE.protein.fa.gz"),
        'datatype':'GENCODE_GTF'
    }
    dc_inputs['RefSeq'] = {
        'gtf': os.path.join(path_of_precisionprodb, "examples", "RefSeq", "RefSeq.gtf.gz"),
        'genome': os.path.join(path_of_precisionprodb, "examples", "RefSeq", "RefSeq.genome.fa.gz"),
        'protein': os.path.join(path_of_precisionprodb, "examples", "RefSeq", "RefSeq.protein.fa.gz"),
        'datatype':'RefSeq'
    }
    dc_inputs['TransDecoder'] = {
        'gtf': os.path.join(path_of_precisionprodb, "examples", "TransDecoder", "TransDecoder.transcripts.fa.transdecoder.genome.gff3.gz"),
        'genome': os.path.join(path_of_precisionprodb, "examples", "TransDecoder", "TransDecoder.genome.fa.gz"),
        'protein': os.path.join(path_of_precisionprodb, "examples", "TransDecoder", "TransDecoder.transcripts.fa.transdecoder.pep.gz",),
        'datatype':'gtf'
    }
    dc_inputs['UniProt'] = {
        'gtf': os.path.join(path_of_precisionprodb, "examples", "Ensembl", "Ensembl.gtf.gz"),
        'genome': os.path.join(path_of_precisionprodb, "examples", "Ensembl", "Ensembl.genome.fa.gz"),
        'protein': os.path.join(path_of_precisionprodb, "examples", "Ensembl", "Ensembl.protein.fa.gz"),
        'UniProt': os.path.join(path_of_precisionprodb, "examples", "UniProt", "UniProt.protein.fa.gz"),
        'datatype':'Ensembl_GTF'
    }

    dc_sqlite = {
        'no_sqlite': '',
        'sqlite_one_step': '',
        'sqlite_two_step': ''
    }
    
    return main_test(dc_variant, dc_inputs, dc_sqlite, output_test, path_of_precisionprodb)


if __name__ == '__main__':
    raise SystemExit(main())
