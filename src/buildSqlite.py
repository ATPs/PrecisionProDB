import pandas as pd
from Bio import SeqIO
import gzip
import os
import time
import pickle
from perChrom import PerChrom
import shutil
import re
import sys
from PrecisionProDB_core import PerGeno
import sqlite3

TEST = False
# TEST = True
if TEST:
    file_genome = '/XCLabServer002_fastIO/examples/GENCODE/GENCODE.genome.fa.gz'
    file_gtf = '/XCLabServer002_fastIO/examples/GENCODE/GENCODE.gtf.gz'
    file_mutations = '/XCLabServer002_fastIO/examples/gnomAD.variant.txt.gz'
    file_protein = '/XCLabServer002_fastIO/examples/GENCODE/GENCODE.protein.fa.gz'
    threads = 10
    outprefix = '/XCLabServer002_fastIO/examples/GENCODE/GENCODE.tsv'
    datatype = 'GENCODE_GTF'
    protein_keyword = 'auto'
    keep_all = True
    file_sqlite = '/XCLabServer002_fastIO/examples/GENCODE/GENCODE.tsv.sqlite'

def insert_chromosomes_using(con, chromosomes_using, chromosomes_genome_description):
    '''
    Store the list of chromosomes in the table "chromosomes_using".
    If the table does not exist, it will be created.
    If a chromosome already exists in the table, update its genome description.
    '''
    cur = con.cursor()
    # Create the table "chromosomes_using" if it doesn't exist
    cur.execute('''
        CREATE TABLE IF NOT EXISTS chromosomes_using (
            id INTEGER PRIMARY KEY AUTOINCREMENT,  -- Add an auto-incrementing primary key
            chromosome TEXT UNIQUE,                  -- Column to store chromosome names with unique constraint
            chromosomes_genome_description TEXT
        )
    ''')
    
    # Upsert (insert or update) the list of chromosomes into the table using executemany for better performance
    cur.executemany('''
        INSERT INTO chromosomes_using (chromosome, chromosomes_genome_description)
        VALUES (?, ?)
        ON CONFLICT(chromosome) DO UPDATE SET
            chromosomes_genome_description=excluded.chromosomes_genome_description
    ''', zip(chromosomes_using, chromosomes_genome_description))

    # Commit the changes
    con.commit()

def read_chromosomes_and_descriptions(cur):
    '''
    Return both chromosomes and their genome descriptions stored in table chromosomes_using.
    The chromosomes and descriptions are returned in two separate lists, maintaining the matching order.
    '''
    try:
        # Query to read the chromosomes and their descriptions from the table
        cur.execute('SELECT chromosome, chromosomes_genome_description FROM chromosomes_using')

        # Fetch all the results
        results = cur.fetchall()
        chromosomes = [result[0] for result in results]
        descriptions = [result[1] for result in results]
        return chromosomes, descriptions
    except:
        print(f"An error occurred while reading from the database")
        return [], []

def insert_df_protein_description(con, df_protein_description, force=False):
    '''
    Store df_protein_description in table protein_description.
    The DataFrame contains a column named "protein_id" that will be used as the unique index for the table.
    If a duplicate protein_id is encountered, keep the newer record.
    If force is True, always create a new table.
    '''
    cur = con.cursor()
    try:
        if force:
            # Drop the existing table if force is True
            cur.execute("DROP TABLE IF EXISTS protein_description")
        
        # Create the table with appropriate column types
        column_definitions = []
        for col in df_protein_description.columns:
            if col == 'protein_id':
                continue
            if pd.api.types.is_integer_dtype(df_protein_description[col].dtype):
                column_definitions.append(f'{col} INTEGER')
            else:
                column_definitions.append(f'{col} TEXT')
        columns = ', '.join(column_definitions)
        
        cur.execute(f'''
            CREATE TABLE IF NOT EXISTS protein_description (
                protein_id TEXT PRIMARY KEY,
                {columns}
            )
        ''')
        
        # Insert or replace rows to ensure newer records overwrite older ones
        cur.executemany(f'''
            INSERT INTO protein_description ({', '.join(df_protein_description.columns)})
            VALUES ({', '.join(['?' for _ in df_protein_description.columns])})
            ON CONFLICT(protein_id) DO UPDATE SET
            {', '.join([f'{col} = excluded.{col}' for col in df_protein_description.columns if col != 'protein_id'])}
        ''', df_protein_description.values.tolist())
        
        # Commit the changes
        con.commit()
    except sqlite3.Error as e:
        print(f"An error occurred while inserting into the database protein_description: {e}")
    finally:
        cur.close()

def get_protein_description(con, protein_ids):
    '''
    Get values from table protein_description based on a list of protein_ids or a single protein_id.
    If protein_ids is a string, return a single record as a DataFrame with one row.
    If protein_ids is a list, return a DataFrame similar to df_protein_description.
    '''
    cur = con.cursor()
    try:
        if isinstance(protein_ids, str):
            # Handle the case where a single protein_id is provided
            query = 'SELECT * FROM protein_description WHERE protein_id = ?'
            cur.execute(query, (protein_ids,))
            row = cur.fetchone()
            if row:
                column_names = [description[0] for description in cur.description]
                df_result = pd.DataFrame([row], columns=column_names)
            else:
                df_result = pd.DataFrame()
        else:
            # Handle the case where a list of protein_ids is provided
            placeholders = ', '.join(['?' for _ in protein_ids])
            query = f'SELECT * FROM protein_description WHERE protein_id IN ({placeholders})'
            cur.execute(query, protein_ids)

            # Fetch all the results and create a DataFrame
            rows = cur.fetchall()
            column_names = [description[0] for description in cur.description]
            df_result = pd.DataFrame(rows, columns=column_names)
        return df_result
    except sqlite3.Error as e:
        print(f"An error occurred while fetching data from the database: {e}")
        return pd.DataFrame()
    finally:
        cur.close()

def insert_df_CDSloc(con, df_CDSloc, force=False):
    '''
    Store df_CDSloc in table CDSloc.
    The DataFrame contains a column named "protein_id" that will be used as an indexed column.
    This column allows duplicate values so that multiple records can be associated with the same protein_id.
    If force is True, always create a new table.
    '''
    cur = con.cursor()
    try:
        if force:
            # Drop the existing table if force is True
            cur.execute("DROP TABLE IF EXISTS CDSloc")
        
        # Create the table with appropriate column types
        column_definitions = []
        for col in df_CDSloc.columns:
            if pd.api.types.is_integer_dtype(df_CDSloc[col].dtype):
                column_definitions.append(f'{col} INTEGER')
            else:
                column_definitions.append(f'{col} TEXT')
        columns = ', '.join(column_definitions)
        
        cur.execute(f'''
            CREATE TABLE IF NOT EXISTS CDSloc (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                {columns},
                FOREIGN KEY (protein_id) REFERENCES protein_description(protein_id)
            )
        ''')
        
        # Create an index on the protein_id column to allow fast lookups
        cur.execute('CREATE INDEX IF NOT EXISTS idx_protein_id ON CDSloc (protein_id)')
        
        # Insert rows into the table
        cur.executemany(f'''
            INSERT INTO CDSloc ({', '.join(df_CDSloc.columns)})
            VALUES ({', '.join(['?' for _ in df_CDSloc.columns])})
        ''', df_CDSloc.values.tolist())
        
        # Commit the changes
        con.commit()
    except sqlite3.Error as e:
        print(f"An error occurred while inserting into the database CDSloc: {e}")
    finally:
        cur.close()

def insert_df_genomicLocs(con, df_genomicLocs, chromosome, force=False):
    '''
    Store df_genomicLocs in table "genomicLocs_" + chromosome.
    The DataFrame contains columns for genomic location information.
    If force is True, always create a new table.
    Additionally, create indexes for genomicLocs_start and genomicLocs_end columns for efficient range queries.
    '''
    table_name = f'genomicLocs_{chromosome}'
    cur = con.cursor()
    try:
        if force:
            # Drop the existing table if force is True
            cur.execute(f"""DROP TABLE IF EXISTS '{table_name}' """)
        
        # Create the table with appropriate column types
        column_definitions = []
        for col in df_genomicLocs.columns:
            if pd.api.types.is_integer_dtype(df_genomicLocs[col].dtype):
                column_definitions.append(f'{col} INTEGER')
            else:
                column_definitions.append(f'{col} TEXT')
        columns = ', '.join(column_definitions)
        
        cur.execute(f'''
            CREATE TABLE IF NOT EXISTS "{table_name}" (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                {columns}
            )
        ''')
        
        # Insert rows into the table
        cur.executemany(f'''
            INSERT INTO "{table_name}" ({', '.join(df_genomicLocs.columns)})
            VALUES ({', '.join(['?' for _ in df_genomicLocs.columns])})
        ''', df_genomicLocs.values.tolist())
        
        # Create indexes for genomicLocs_start and genomicLocs_end to enable efficient range queries
        cur.execute(f'CREATE INDEX IF NOT EXISTS idx_genomicLocs_start_end ON "{table_name}" (genomicLocs_start, genomicLocs_end)')
        
        # Commit the changes
        con.commit()
    except sqlite3.Error as e:
        print(f"An error occurred while inserting into the database '{table_name}': {e}")
    finally:
        cur.close()

def get_CDSloc(con, protein_ids):
    '''
    Get values from table CDSloc based on a list of protein_ids or a single protein_id.
    If protein_ids is a string, return all records for that protein_id as a DataFrame.
    If protein_ids is a list, return a DataFrame containing all records matching any of the protein_ids.
    '''
    cur = con.cursor()
    try:
        if isinstance(protein_ids, str):
            # Handle the case where a single protein_id is provided
            query = 'SELECT * FROM CDSloc WHERE protein_id = ?'
            cur.execute(query, (protein_ids,))
            rows = cur.fetchall()
            if rows:
                column_names = [description[0] for description in cur.description]
                df_result = pd.DataFrame(rows, columns=column_names)
            else:
                df_result = pd.DataFrame()
        else:
            # Handle the case where a list of protein_ids is provided
            placeholders = ', '.join(['?' for _ in protein_ids])
            query = f'SELECT * FROM CDSloc WHERE protein_id IN ({placeholders})'
            cur.execute(query, protein_ids)

            # Fetch all the results and create a DataFrame
            rows = cur.fetchall()
            column_names = [description[0] for description in cur.description]
            df_result = pd.DataFrame(rows, columns=column_names)
        
        if not df_result.empty:
            df_result = df_result.sort_values(by = 'id')
            # drop column "id"
            df_result = df_result.drop(columns = ['id'])
        return df_result
    except sqlite3.Error as e:
        print(f"An error occurred while fetching data from the database: {e}")
        return pd.DataFrame()
    finally:
        cur.close()

def get_genomicLocs(con, protein_ids, chromosome=None):
    '''
    Get values from the appropriate genomicLocs table(s) based on a list of protein_ids or a single protein_id.
    If chromosome is provided, search the table specific to that chromosome.
    If chromosome is None, search all genomicLocs tables.
    '''
    cur = con.cursor()
    try:
        df_result = pd.DataFrame()
        if chromosome:
            # Search a specific chromosome table
            table_name = f'genomicLocs_{chromosome}'
            if isinstance(protein_ids, str):
                query = f'SELECT * FROM "{table_name}" WHERE protein_id = ?'
                cur.execute(query, (protein_ids,))
                rows = cur.fetchall()
            else:
                placeholders = ', '.join(['?' for _ in protein_ids])
                query = f'SELECT * FROM "{table_name}" WHERE protein_id IN ({placeholders})'
                cur.execute(query, protein_ids)
                rows = cur.fetchall()
            if rows:
                column_names = [description[0] for description in cur.description]
                df_result = pd.DataFrame(rows, columns=column_names)
        else:
            # Search across all genomicLocs tables
            cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'genomicLocs_%'")
            tables = cur.fetchall()
            for table in tables:
                table_name = table[0]
                if isinstance(protein_ids, str):
                    query = f'SELECT * FROM "{table_name}" WHERE protein_id = ?'
                    cur.execute(query, (protein_ids,))
                    rows = cur.fetchall()
                else:
                    placeholders = ', '.join(['?' for _ in protein_ids])
                    query = f'SELECT * FROM "{table_name}" WHERE protein_id IN ({placeholders})'
                    cur.execute(query, protein_ids)
                    rows = cur.fetchall()
                if rows:
                    column_names = [description[0] for description in cur.description]
                    df_table = pd.DataFrame(rows, columns=column_names)
                    df_result = pd.concat([df_result, df_table], ignore_index=True)
        
        if not df_result.empty:
            df_result = df_result.sort_values(by='id').reset_index(drop=True)
            # drop column "id"
            df_result = df_result.drop(columns=['id'])
        return df_result
    except sqlite3.Error as e:
        print(f"An error occurred while fetching data from the database: {e}")
        return pd.DataFrame()
    finally:
        cur.close()

def get_protein_id_from_genomicLocs(con, chromosome, pos):
    '''
    Given a chromosome and a position, return the protein_ids from the genomicLocs table of that chromosome,
    where genomicLocs_start <= pos <= genomicLocs_end.
    '''
    table_name = f'genomicLocs_{chromosome}'
    pos = int(pos)
    cur = con.cursor()
    try:
        # Query to find rows where the position falls between genomicLocs_start and genomicLocs_end
        query = f'SELECT protein_id FROM "{table_name}" WHERE genomicLocs_start <= ? AND genomicLocs_end >= ?'
        cur.execute(query, (pos, pos))
        rows = cur.fetchall()
        # Extract the protein_ids from the result rows
        protein_ids = [row[0] for row in rows]
        return protein_ids
    except sqlite3.Error as e:
        print(f"An error occurred while fetching data from the database: {e}")
        return []
    finally:
        cur.close()

def create_sqlite(file_sqlite, file_genome, file_gtf, file_protein, outprefix, datatype, protein_keyword, threads=None, keep_all=False):
    '''
    Create a SQLite database containing genomic and protein information.

    Parameters:
    - file_sqlite: Path to the SQLite output file.
    - file_genome: Path to the genome fasta file.
    - file_gtf: Path to the GTF file with gene annotations.
    - file_protein: Path to the protein sequences file.
    - outprefix: Prefix for output files generated.
    - datatype: Type of input data (e.g., GENCODE_GTF, GENCODE_GFF3).
    - protein_keyword: Keyword for filtering protein sequences.
    - keep_all: Boolean to indicate if all data should be kept.
    '''

    pergeno = PerGeno(
        file_genome = file_genome, 
        file_gtf=file_gtf, 
        file_mutations = '', 
        file_protein=file_protein, 
        threads=threads, 
        outprefix=outprefix, 
        datatype=datatype, 
        protein_keyword=protein_keyword, 
        keep_all = keep_all
        )
    # print(pergeno.__dict__)
    pergeno.splitInputByChromosomes()
    
    # process each chromosome
    chromosomes = pergeno.chromosomes # chromosomes with proteins
    chromosomes_genome = pergeno.chromosomes_genome # all chromosomes in genome
    chromosomes_genome_description = pergeno.chromosomes_genome_description

    # create a sqlite3 file
    con = sqlite3.connect(file_sqlite)
    cur = con.cursor()

    # create a table to store chromosomes
    insert_chromosomes_using(con, chromosomes_using=chromosomes_genome, chromosomes_genome_description=chromosomes_genome_description)

    # run get df_transcript2
    for chromosome in chromosomes:
        fchr_genome = os.path.join(pergeno.tempfolder, chromosome + '.genome.fa')
        fchr_gtf = os.path.join(pergeno.tempfolder, chromosome + '.gtf')
        fchr_protein = os.path.join(pergeno.tempfolder, chromosome + '.proteins.fa')
        perchrom = PerChrom(file_genome = fchr_genome,
            file_gtf = fchr_gtf,
            file_mutations = '',
            file_protein = fchr_protein,
            threads = 1,
            outprefix = outprefix,
            datatype = pergeno.datatype,
            chromosome = chromosome)
        
        df_transcript2 = perchrom.df_transcript2
        # store protein_description, AA_seq to table protein_description
        df_protein_description = df_transcript2[[i for i in df_transcript2.columns if i not in ['CDSloc']]].copy()
        df_protein_description['genomicLocs'] = df_protein_description['genomicLocs'].astype(str)
        df_protein_description = df_protein_description.reset_index()
        insert_df_protein_description(con, df_protein_description)

        # store CDSloc (location of CDS) to table CDSloc
        df_CDSloc = df_transcript2[['CDSloc']].copy()
        df_CDSloc = df_CDSloc.explode('CDSloc')
        df_CDSloc['CDSloc_start'] = df_CDSloc['CDSloc'].str[0] + 1 # 1-based
        df_CDSloc['CDSloc_end'] = df_CDSloc['CDSloc'].str[1]
        df_CDSloc = df_CDSloc.reset_index()
        df_CDSloc['seqname'] = chromosome
        del df_CDSloc['CDSloc']
        insert_df_CDSloc(con, df_CDSloc, force=False)

        # store genomicLocs (location of CDSplus) to based on each chromosome. genomicLocs is the loc of CDSplus, not CDS
        df_genomicLocs = df_transcript2[['genomicLocs']].copy()
        df_genomicLocs = df_genomicLocs.explode('genomicLocs')
        df_genomicLocs['genomicLocs_start'] = df_genomicLocs['genomicLocs'].str[0] + 1 # 1-based
        df_genomicLocs['genomicLocs_end'] = df_genomicLocs['genomicLocs'].str[1]
        df_genomicLocs = df_genomicLocs.reset_index()
        df_genomicLocs['seqname'] = chromosome
        del df_genomicLocs['genomicLocs']
        insert_df_genomicLocs(con, df_genomicLocs, chromosome, force=False)

    print('building sqlite done')
    con.close()

def get_connection(file_sqlite):
    '''Get a connection to the SQLite database.'''
    try:
        conn = sqlite3.connect(file_sqlite)
        print(f'Connected to SQLite database: {file_sqlite}')
        return conn
    except sqlite3.Error as e:
        print(f'Error connecting to SQLite database: {e}')
        sys.exit(1)

description = '''create a sqlite file for PrecisionProDB.'''

if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-q', '--sqlite', help='SQLite file with CDS and exon annotations', required=True)
    parser.add_argument('-g', '--genome', help='Path to the genome fasta file', required=True)
    parser.add_argument('-p', '--protein', help='Path to the protein sequences file', required=True)
    parser.add_argument('-f', '--gtf', help='Path to the GTF file with gene annotations', required=True)
    parser.add_argument('-o', '--out', help='Output prefix. Two files will be output: one with annotation for mutated transcripts and one with the protein sequences. {out}.aa_mutations.csv, {out}.mutated_protein.fa', default=None)
    parser.add_argument('-a', '--datatype', help='''input datatype, could be GENCODE_GTF, GENCODE_GFF3, RefSeq, Ensembl_GTF or gtf. default "gtf". Ensembl_GFF3 is not supported. ''', default='gtf', type=str, choices=['GENCODE_GTF', 'GENCODE_GFF3','RefSeq','Ensembl_GTF','gtf'])
    parser.add_argument('-k', '--protein_keyword', help='Keyword for filtering protein sequences', default='auto')
    parser.add_argument('--keep_all', help='Keep all data. Default is False.', type=bool, default=False)
    
    if TEST:
        args = parser.parse_args(['-q', file_sqlite, '-g', file_genome, '-p', file_protein, '-f', file_gtf, '-o', outprefix, '-a', 'GENCODE_GTF', '-k', protein_keyword, '--keep_all', str(keep_all)])

    else:
        args = parser.parse_args()
    
    create_sqlite(
        file_sqlite=args.sqlite,
        file_genome=args.genome,
        file_gtf=args.gtf,
        file_protein=args.protein,
        outprefix=args.out,
        datatype=args.datatype,
        protein_keyword=args.protein_keyword,
        keep_all=args.keep_all,
    )
    
    print('SQLite creation complete.')
