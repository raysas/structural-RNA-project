from setuptools import setup, find_packages
import os
import re

# Read version from __init__.py without importing
def get_version():
    init_path = os.path.join(os.path.dirname(__file__), 'src', 'rna_score', '__init__.py')
    with open(init_path, 'r') as f:
        content = f.read()
        version_match = re.search(r"^__version__\s*=\s*['\"]([^'\"]*)['\"]", content, re.M)
        if version_match:
            return version_match.group(1)
        raise RuntimeError("Unable to find version string.")
    
# def get_requirements():
#     req_path = os.path.join(os.path.dirname(__file__), 'requirements.txt')
#     with open(req_path, 'r') as f:
#         requirements = []
#         for line in f:
#             line = line.strip()
#             if line and not line.startswith('#'):
#                 requirements.append(line)
#         return requirements

__version__ = get_version()

setup(
    name='rna_score',
    version=__version__,
    description='RNA structure scoring library',
    author='Yazid Hoblos, Joelle Assy, Denys Buryi, Raul Duran De Alba, Rayane Adam',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    py_modules=['access_rna_structures'],
    scripts=['src/access_rna_structures.py'],
    install_requires=['numpy'],
    entry_points={
        'console_scripts': [
            'rna-score = rna_score.cli:main',
        ],
    },
)
