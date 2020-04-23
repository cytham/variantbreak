from setuptools import setup, find_packages
from variantbreak.version import __version__
import os

current_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(current_dir, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='variantbreak',
    version=__version__,
    packages=find_packages(),
    scripts=['variantbreak/variantbreak'],
    url='https://github.com/cytham/variantbreak',
    download_url='https://github.com/cytham/variantbreak/releases',
    license='gpl-3.0',
    author='Tham Cheng Yong',
    author_email='chengyong.tham@u.nus.edu',
    description='Structural variant analyzer for data visualization on VariantMap',
    keywords=['variantbreak', 'structural variant', 'sv', 'breakend', 'annotation', 'nanovar', 'nanopore', 'long read',
              'variantmap, filter'],
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=['numpy>=1.17.3', 'scipy>=1.2.1', 'biopython>=1.74', 'pybedtools>=0.8.0', 'matplotlib>=2.2.3',
                      'tensorflow>=2.0.0', 'natsort>=6.2.0', 'progress>=1.4', 'pysam>=0.15.3'],
    python_requires='>=3.6',
    classifiers=[
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
