"""
Main annotation function of VariantBreak.

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

from __future__ import absolute_import

import os
import sys
import logging
from datetime import datetime

import csv
import pandas as pd
from .valid import valid
from .parse import variant_parse, annotation_parse, filter_parse


def variantbreak(
        variant_file_dir,
        annotation_file_dir=None,
        filter_file_dir=None,
        # ignore_duplicate=False,
        del_annotate_size=20000,
        merge_buffer=50,
        promoter_size=1000,
        max_annotation=3,
        cluster_sample=False,
        auto_filter=False,
        quiet=False
):
    """ Returns a dataframe with metadata

- variant_file_dir (str; required): Path to a single variant file or directory containing
    variant files in VCF or BED formats. VCF files needs to be generated by NanoVar long-read
    SV caller (v1.3.5 or above). BED files need to be tab-delimited with the essential three
    columns: chrom, chromStart and chromEnd (column 4 will be ignored).
- annotation_file_dir (str; optional): Path to a single annotation file or directory containing
    annotation files of GTF/GFF or BED formats.
- filter_file_dir (str; optional): Path to a single filter file or directory containing filter
    files of BED format.
- del_annotate_size (int; default 20000): Deletions with sizes lower or equal to this value will
    have the entire deleted region annotated. Any genes that intersect with the deleted region
    will be included as annotation. On the contrary, if deletion size is greater than this value,
    only the deletion breakends will be annotated, omitting any annotation of the middle deleted
    region. In other words, only genes intersecting with the deletion breakends will be included
    as annotation. This is to reduce false annotations due to false large deletions. Note that
    the value to be set is an absolute deletion size, do not use minus '-'. Use value '-1' to
    disable this threshold and annotate all deleted regions despite of size.
- merge_buffer (int; default 50): Nucleotide length buffer for SV breakend clustering
- promoter_size (int; default 1000): Length in base-pairs upstream of TSS to define promoter region
- max_annotation (int; default 3): Maximum number of annotation entries to be recorded in the
    dataframe for each SV
- cluster_sample (bool; default False): If True, performs hierarchical clustering on samples.
- auto_filter (bool; default False): If True, automatically removes variants that intersected
    with all filter BED files.
- quiet (bool; default False): Hide verbose if True.


Usage example:

# Import variantbreak function from variantbreak package
from variantbreak import variantbreak

# Run variantbreak on your samples in sample_dir, together with annotation and filter files, and output dataframe to df
df = variantbreak("/path/to/sample_dir/",
                  "/path/to/annotation_dir/",
                  "/path/to/filter_dir/")

# To save data to files
# Import write_to_file from variantbreak package
from variantbreak import write_to_files

# Specify dataframe variable, output file path and prefix, and delimiter
write_to_files(df,
               "/path/to/output_prefix",
               sep="\t")

# Output files
output.h5 - HDF5 file required for data visualization by VariantMap
output.csv - CSV file for data viewing, separated by the delimiter set by user
legend.txt - File containing the legend of the sample labels used in the dataframe

    """

    # Set verbose
    if quiet:
        sys.stdout = open(os.devnull, 'w')

    # Disable del_annotate_size if it is -1
    if int(del_annotate_size) == -1:
        del_annotate_size = float('inf')

    # Check file integrity and right formatting for input files
    file_dict = valid(variant_file_dir, annotation_file_dir, filter_file_dir)

    # Print and record time started
    now = datetime.now()
    now_str = now.strftime("[%d/%m/%Y %H:%M:%S]")
    print(now_str, "- VariantBreak started")
    logging.info('VariantBreak started')

    # Parse annotation files
    annotation_dict, gene_info_dict, label_col, annote_name_dict = annotation_parse(file_dict, promoter_size)

    # Parse filter files
    filter_dict, filter_col, filter_name_dict = filter_parse(file_dict)

    # Parse and intersect variant files and generate master dataframe
    df, sample_name_dict = variant_parse(file_dict,
                                         del_annotate_size,
                                         merge_buffer,
                                         annotation_dict,
                                         gene_info_dict,
                                         label_col,
                                         filter_dict,
                                         filter_col,
                                         max_annotation,
                                         cluster_sample,
                                         auto_filter
                                         )

    # Gather sample labels
    sample_names = [x for x in sample_name_dict.keys()]

    # Add metadata to dataframe
    df.metadata = ''  # To avoid pandas UserWarning
    df.metadata = {'sample_names': sample_names, 'sample': sample_name_dict, 'annotation': annote_name_dict,
                   'filter': filter_name_dict}

    # Print and record time finished
    now = datetime.now()
    now_str = now.strftime("[%d/%m/%Y %H:%M:%S]")
    print(now_str, "- VariantBreak ended")
    logging.info('VariantBreak ended')
    return df


def write_to_files(df, file_path, sep=','):

    # Save dataframe and metadata to HDF5
    filename = file_path  # , _ = os.path.splitext(file_path)
    store = pd.HDFStore(filename + '.h5', mode='w')
    store.put('dataset', df)
    metadata = df.metadata
    store.get_storer('dataset').attrs.metadata = metadata
    store.close()

    # Save dataframe to CSV
    df2 = df.copy()
    for name in metadata['sample_names']:
        df2[name].replace({0.0: "-", 0.214: "DEL", 0.357: "INV", 0.500: "INS", 0.643: "BND", 0.786: "DUP", 0.929: "UKN"},
                          inplace=True)
    df2.to_csv(filename + '.csv', index=True, sep=sep, index_label='Index')

    # Save label names to file_path in a legend text file for variant and filter files
    w = open(os.path.join(os.path.dirname(filename), 'legend.txt'), 'w')
    w.write('#Label' + sep + 'File_path\n')
    write_name = csv.writer(w, delimiter=sep)
    legend = df.metadata['sample']
    legend.update(df.metadata['filter'])
    legend.update(df.metadata['annotation'])
    for key, val in legend.items():
        write_name.writerow([key, val])
    w.close()
