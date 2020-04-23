"""
Functions for parsing and verifying input parameters.

Copyright (C) 2020 Tham Cheng Yong

This file is part of VariantBreak.

VariantBreak is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

VariantBreak is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with VariantBreak.  If not, see <https://www.gnu.org/licenses/>.
"""

import sys
import argparse

import pandas as pd
from variantbreak import __version__


# Parse input
def input_parser(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="VariantBreak is a structural variant (SV) breakend analyzer that \
intersects and merges \nSV breakends from NanoVar VCF files or variant BED files for visualization on VariantMap \nor \
summarized into a CSV file. It also annotates and filters all SVs across all samples \naccording to user input GTF/GFF/BED \
files.",
                                     formatter_class=argparse.RawTextHelpFormatter)

    def restrict_int(f):
        try:
            if int(f) < -1:
                raise argparse.ArgumentTypeError("%r needs to be positive or -1 to disable" % (int(f),))
            else:
                return f
        except ValueError:
            raise argparse.ArgumentTypeError("%r needs to be an integer" % (f,))

    def one_character(s):
        df = pd.DataFrame()
        _sep = bytes(s, "utf-8").decode("unicode_escape")
        try:
            df.to_csv(sep=_sep)
            return _sep
        except TypeError:
            raise argparse.ArgumentTypeError("%r needs to be a single character" % (s,))

    parser.add_argument("v_file_dir", type=str,
                        metavar="[variant_path]",
                        help="""path to single variant file or directory containing variant files 
of VCF (NanoVar-v1.3.6 or above) or BED formats. Formats: .vcf/.bed""")

    parser.add_argument("wk_dir", type=str,
                        metavar="[working_directory]",
                        help="""path to working directory. Directory will be created if it does not 
exist.""")

    parser.add_argument("-a", "--annotation_file_dir", type=str,
                        metavar="path",
                        help="""path to single annotation file or directory containing annotation 
files of GTF/GFF or BED formats. Formats: .gtf/.gff/.gff3/.bed""")

    parser.add_argument("-f", "--filter_file_dir", type=str,
                        metavar="path",
                        help="""path to single filter file or directory containing filter files of 
BED format. Format: .bed""")

    parser.add_argument("-d", "--del_annotate_size", type=restrict_int, metavar="int",
                        default=20000,
                        help="""Deletions with sizes lower or equal to this value will have the 
entire deleted region annotated. Any genes that intersect with 
the deleted region will be included as annotation. On the contrary, 
if deletion size is greater than this value, only the deletion 
breakends will be annotated, omitting any annotation of the middle 
deleted region. In other words, only genes intersecting with the 
deletion breakends will be included as annotation. This is done to 
reduce false annotations due to false large deletions. Note that 
the value to be set is an absolute deletion size, do not use minus 
'-'. Use value '-1' to disable this threshold and annotate all 
deleted regions despite of size. [20000]""")

    parser.add_argument("-b", "--merge_buffer", type=int, metavar="int",
                        default=50,
                        help="""nucleotide length buffer for SV breakend clustering [50]""")

    parser.add_argument("-p", "--promoter_size", type=int, metavar="int",
                        default=1000,
                        help="""length in base-pairs upstream of TSS to define promoter region 
[1000]""")

    parser.add_argument("-m", "--max_annotation", type=int, metavar="int",
                        default=3,
                        help="""maximum number of annotation entries to be recorded in the 
dataframe for each SV [3]]""")

    parser.add_argument("-s", "--sep", type=one_character, metavar="str",
                        default=',',
                        help="""single character field delimiter for output dataframe CSV file 
(e.g. '\\t' for tab or ',' for comma) [,]""")

    parser.add_argument("-n", "--filename", type=str, metavar="str",
                        default="output",
                        help="File name prefix of output files [output]")

    parser.add_argument("-v", "--version", action='version',
                        version=__version__,
                        help="show version and exit")

    parser.add_argument("-q", "--quiet", action='store_true',
                        help="hide verbose")

    args = parser.parse_args(args)
    return args
