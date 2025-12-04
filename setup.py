from setuptools import setup, find_packages

setup(
    name='rna_score',
    version='0.1.0',
    description='RNA structure scoring library',
    author='Your Name',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=['numpy'],
    entry_points={
        'console_scripts': [
            'rna-score = rna_score.cli:main',
        ],
    },
)
