import os
import argparse
import subprocess
import time


def run_command(command):
    """
    Run a system command using subprocess and print the output.
    """
    print(f"Running command: {command}")
    start_time = time.time()
    result = subprocess.run(command, shell=True)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Command finished in {elapsed_time:.2f} seconds\n\n")
    if result.returncode != 0:
        print(f"Error running command: {command}\n\n")
        return command, 'failed'
    else:
        return command, 'completed'

def get_cmd_to_run(key_input, key_variant, sqlite_key):
    '''
    '''
    if key_input == 'UniProt' and key_variant != 'str':
        file_UniProt = dc_inputs[key_input]['UniProt']
        cmd_UniProt = f' -U {file_UniProt} -t 4 -D Uniprot '
    else:
        cmd_UniProt = ' -t 4 '
    
    folder_work = os.path.join(output_test, key_input, key_variant, sqlite_key)
    if not os.path.exists(folder_work):
        os.makedirs(folder_work)
    file_mutation = dc_variant[key_variant]
    file_genome = dc_inputs[key_input]['genome']
    file_protein = dc_inputs[key_input]['protein']
    file_gtf = dc_inputs[key_input]['gtf']
    datatype = dc_inputs[key_input]['datatype']
    file_sqlite = f'{output_test}/{key_input}/{key_input}.sqlite'
    output = os.path.join(folder_work, f'{key_input}.{key_variant}.{sqlite_key}')
    
    print(f'running test with key_input: {key_input}, key_variant: {key_variant}, sqlite_key: {sqlite_key}')
    if sqlite_key == 'no_sqlite':
        print("Running test: without use sqlite file")
        cmd = f'cd {folder_work} &&  python {path_of_precisionprodb}/src/PrecisionProDB.py -m {file_mutation} -g {file_genome} -p {file_protein} -f {file_gtf} -o {output} -a {datatype} --PEFF {cmd_UniProt}'
    elif sqlite_key == 'sqlite_one_step':
        print("Running test: use SQLite file as intermediate file")
        if os.path.exists(file_sqlite):
            print(f"{file_sqlite} exists. Removing it.")
            os.remove(file_sqlite)
        cmd = f'cd {folder_work} &&  python {path_of_precisionprodb}/src/PrecisionProDB.py -m {file_mutation} -g {file_genome} -p {file_protein} -f {file_gtf} -o {output} -a {datatype} --PEFF {cmd_UniProt} -S {file_sqlite}'
    elif sqlite_key == 'sqlite_two_step':
        print("Running test: Generate SQLite file in advance and use SQLite")
        if os.path.exists(file_sqlite):
            print(f"{file_sqlite} exists. Removing it.")
            os.remove(file_sqlite)
        cmd = f'cd {folder_work} && python {path_of_precisionprodb}/src/buildSqlite.py -S {file_sqlite} -g {file_genome} -p {file_protein} -f {file_gtf}  -a {datatype} && python {path_of_precisionprodb}/src/PrecisionProDB.py -m {file_mutation}  -o {output} -a {datatype} --PEFF {cmd_UniProt} -S {file_sqlite}'
    
    return cmd

def main(dc_variant, dc_inputs, dc_sqlite, output_test):

    ls_cmd = []
    ls_results = []
    for key_input in dc_inputs:
        for key_variant in dc_variant:
            for sqlite_key in dc_sqlite:
                if sqlite_key == 'no_sqlite' and key_variant == 'str':
                    continue
                cmd = get_cmd_to_run(key_input, key_variant, sqlite_key)
                ls_cmd.append(cmd)
                ls_results.append(run_command(cmd))

    # save ls_results to file
    file_job_summary = os.path.join(output_test, 'test_running_summary.txt')
    with open(file_job_summary, 'w') as f:
        for i, result in enumerate(ls_results):
            f.write(f'{ls_cmd[i]}\n{result}\n\n')


if __name__ == '__main__':
    description = """
    Run PrecisionProDB tests.
    Usage: test.py -s PATH_OF_PRECISIONPRODB [-o OUTPUT_TEST]
    Make sure Python is in your PATH

    Options:
        -s PATH_OF_PRECISIONPRODB  Path to PrecisionProDB. Will use PATH_OF_PRECISIONPRODB/src/ to find the scripts
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
    path_of_precisionprodb = (path_of_precisionprodb if path_of_precisionprodb else os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    output_test = (
        args.output if args.output else os.path.join(path_of_precisionprodb, "test_output")
    )

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
    
    main(dc_variant, dc_inputs, dc_sqlite, output_test)







