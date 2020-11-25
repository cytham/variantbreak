## VariantBreak - Structural variant analyzer for data visualization on VariantMap
[![Build Status](https://travis-ci.com/cytham/variantbreak.svg?branch=master)](https://travis-ci.com/cytham/variantbreak)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/variantbreak)](https://pypi.org/project/variantbreak/)
[![PyPI versions](https://img.shields.io/pypi/v/variantbreak)](https://pypi.org/project/variantbreak/)
[![Conda](https://img.shields.io/conda/v/bioconda/variantbreak)](https://anaconda.org/bioconda/variantbreak)
[![Github release](https://img.shields.io/github/v/release/cytham/variantbreak?include_prereleases)](../../releases)
[![PyPI license](https://img.shields.io/pypi/l/variantbreak)](./LICENSE.txt)

VariantBreak is a python package that integrates all structural variants (SVs) from a cohort of 
[NanoVar](https://github.com/cytham/nanovar) VCF files or variant BED files for visualization on [VariantMap](https://github
.com/cytham/variantmap) or summarized into a CSV file. It also annotates and filters all SVs across all samples according to
 user input GTF/GFF/BED files. 

### Basic capabilities
* Intersects and merges all SV breakends from a sample cohort using [NanoVar](https://github.com/cytham/nanovar) VCF files 
(NanoVar-v1.3.6 or above) or variant BED files.
* Annotates each SV according to input GTF/GFF files or BED annotation files.
* Filters SVs by adding a "HIT" or "MISS" label according to input BED filter files.
* Creates a master pandas dataframe to store all data. 
* Creates a HDF5 file containing the master dataframe and some metadata which can be graphically visualized on VariantMap
 within Dash Bio.

## Getting Started

### Quick run

##### Command-line usage:
```
variantbreak [Options] -a annotation.gff3 -f filter.bed variant_path working_dir 
```

| Parameter | Argument | Comment |
| :--- | :--- | :--- |
| `-a` | annotation.gff3 | path to single annotation file or directory containing annotation files of GTF/GFF or BED formats |
| `-f` | filter.bed | path to single filter file or directory containing filter files of BED format|
| - | variant_path | path to single variant file or directory containing variant files of VCF or BED formats|
| - | working_dir | path to working directory |

##### Python console usage:
```
# Import variantbreak function from variantbreak package
from variantbreak import variantbreak

# Run variantbreak on your samples with annotation and filter files
df = variantbreak("/path/to/sample_dir/",
                  "/path/to/annotation_dir/",
                  "/path/to/filter_dir/")


# To save data to files
# Import write_to_file from variantbreak package
from variantbreak import write_to_files

# Specify dataframe variable, output file path and prefix, and delimiter of choice
write_to_files(df,
               "/path/to/output_prefix",
               sep="\t")

```
#### Output
| Output file | Comment |
| :--- | :--- |
| output.h5 | HDF5 file required for data visualization by VariantMap |
| output.csv | CSV file for data viewing, separated by the delimiter set by user |
| legend.txt | File containing the legend of the sample labels used in analysis|

For more information, see [wiki](https://github.com/cytham/variantbreak/wiki).

### Operating system: 
* Linux (x86_64 architecture, tested in Ubuntu 16.04)

### Installation:
There are three ways to install VariantBreak:
#### Option 1: Conda (Recommended)
```
# Installing from bioconda automatically installs all dependencies 
conda install -c bioconda variantbreak
```
#### Option 2: Pip (See dependencies below)
```
# Installing from PyPI requires own installation of dependencies, see below
pip install variantbreak
```
#### Option 3: GitHub (See dependencies below)
```
# Installing from GitHub requires own installation of dependencies, see below
git clone https://github.com/cytham/variantbreak.git 
cd variantbreak
pip install .
```

### Installation of dependencies
* bedtools >=2.26.0 (required to be in PATH by pybedtools)
* pybedtools >=0.8.1
* pandas >=1.0.3
* tables >=3.6.1
* fastcluster >=1.1.26

##### 1. _bedtools_
Please visit [here](https://bedtools.readthedocs.io/en/latest/content/installation.html) for instructions to install.

##### 2. _pybedtools_
Please visit [here](https://daler.github.io/pybedtools/main.html) for instructions to install.

##### 3. _pandas_
Please visit [here](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html) for instructions to install.

##### 4. _tables_
```
pip install tables
```
or
```
conda install -c conda-forge pytables
```

##### 5. _fastcluster_
```
pip install fastcluster
```
or
```
conda install -c conda-forge fastcluster
```

## Documentation
See [wiki](https://github.com/cytham/variantbreak/wiki) for more information.

## Versioning
See [CHANGELOG](./CHANGELOG.txt)

## Citation
Not available

## Author

* **Tham Cheng Yong** - [cytham](https://github.com/cytham)

## License

VariantBreak is licensed under GNU General Public License - see [LICENSE.txt](./LICENSE.txt) for details.

## Limitations
* Current version only allows input of VCF files generated by NanoVar. We will create a format adaptor in future versions to
 encompass VCF files generated by other SV callers.
 
* Processing speed of large sample cohorts has not been tested. Currently, it takes about 30 minutes to process about 100,000
 merged SVs. 
