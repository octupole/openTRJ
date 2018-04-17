"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding

# Get the long description from the README file
"""
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
"""
# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='openTRJ',  # Required
    version='0.2',  # Required
    description='the openTRJ python scripts',  # Required
    author='Massimo Marchi',  # Optional
    packages=['Micelles'],
    install_requires=['numpy >= 1.10','ujson >=1.2','keyring','paramiko'],
    entry_points = {
        'console_scripts': [
            'Rg = Rg:main',
            'Voro = Voro:main'
        ]
    }
#    scripts=['bin/Rg.py','bin/Voro.py']
)
