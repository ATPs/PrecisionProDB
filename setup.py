from setuptools import setup, find_packages
import os
folder = os.path.dirname(os.path.abspath(__file__))

# Read version from version file
def get_version():
    version_file = os.path.join(folder, 'src', 'precisionprodb', 'version')
    with open(version_file, 'r') as f:
        return f.read().strip()

setup(
    name="precisionprodb",  # Your package name
    version=get_version(),  # Package version read from package version file
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    package_data={'precisionprodb': ['version']},
    install_requires=[  # Your dependencies
        'biopython>=1.78',
        'pandas>=1.0.5',
        'numpy>=1.18.5',
    ],
    entry_points={  # Define terminal commands for each script
        'console_scripts': [
            'buildSqlite=precisionprodb.buildSqlite:main',
            'buildPeptideSqlite=precisionprodb.peptideSqlite:main',
            'downloadHuman=precisionprodb.downloadHuman:main',
            'extractMutatedUniprot=precisionprodb.extractMutatedUniprot:main',
            'generatePEFFoutput=precisionprodb.generatePEFFoutput:main',
            'perChrom=precisionprodb.perChrom:main',
            'perChromSqlite=precisionprodb.perChromSqlite:main',
            'PrecisionProDB=precisionprodb.PrecisionProDB:main',
            'PrecisionProDB_core=precisionprodb.PrecisionProDB_core:main',
            'PrecisionProDB_Sqlite=precisionprodb.PrecisionProDB_Sqlite:main',
            'PrecisionProDB_test=precisionprodb.PrecisionProDB_test:main',
            'PrecisionProDB_vcf=precisionprodb.PrecisionProDB_vcf:main',
            'vcf2mutation=precisionprodb.vcf2mutation:main',
        ],
    },
    long_description=open(folder +'/README.md').read(),  # Read the README as long description
    long_description_content_type='text/markdown',
    author="Xiaolong Cao",
    author_email="ATPs@outlook.com",
    url="https://github.com/ATPs/PrecisonProDB",
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3.0',  # Adjust to your license
        'Operating System :: OS Independent',
    ],
)
