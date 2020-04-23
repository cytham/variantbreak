"""
Functions for validating input files.

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

import os
from distutils.version import StrictVersion


def valid(
        variant_file_dir,
        annotation_file_dir,
        filter_file_dir,
):
    """Return a dictionary with file path or name list from different input directories.

Keyword arguments:

- variant_file_dir (path; required): Variant file or path to directory containing all variant
    files of .vcf or .bed formats.
- annotation_file_dir (path): Annotation file or path to directory containing all annotation
    files of .gff, .gff3, .gtf, .bed formats.
- filter_file_dir (path): Filter file or path to directory containing all filter files of .bed format.
    """

    # Get file lists
    if os.path.isfile(variant_file_dir):
        file_dict = {'variant_list': [variant_file_dir],
                     'variant_name_list': [os.path.basename(variant_file_dir)]}
    else:
        file_dict = {'variant_list': sorted([os.path.join(variant_file_dir, file)
                                             for file in os.listdir(variant_file_dir)
                                             if os.path.isfile(os.path.join(variant_file_dir, file))]),
                     'variant_name_list': sorted([file for file in os.listdir(variant_file_dir)
                                                  if os.path.isfile(os.path.join(variant_file_dir, file))])}

    if annotation_file_dir:
        if os.path.isfile(annotation_file_dir):
            file_dict['label_list'] = [annotation_file_dir]
            file_dict['label_name_list'] = [os.path.basename(annotation_file_dir)]
        else:
            file_dict['label_list'] = sorted([os.path.join(annotation_file_dir, file)
                                              for file in os.listdir(annotation_file_dir)
                                              if os.path.isfile(os.path.join(annotation_file_dir, file))])
            file_dict['label_name_list'] = sorted([file for file in os.listdir(annotation_file_dir)
                                                   if os.path.isfile(os.path.join(annotation_file_dir, file))])
    else:
        file_dict['label_list'], file_dict['label_name_list'] = [], []

    if filter_file_dir:
        if os.path.isfile(filter_file_dir):
            file_dict['filter_list'] = [filter_file_dir]
            file_dict['filter_name_list'] = [os.path.basename(filter_file_dir)]
        else:
            file_dict['filter_list'] = sorted([os.path.join(filter_file_dir, file)
                                               for file in os.listdir(filter_file_dir)
                                               if os.path.isfile(os.path.join(filter_file_dir, file))])
            file_dict['filter_name_list'] = sorted([file for file in os.listdir(filter_file_dir)
                                                    if os.path.isfile(os.path.join(filter_file_dir, file))])
    else:
        file_dict['filter_list'], file_dict['filter_name_list'] = [], []

    # Validate files
    for file in file_dict['variant_list']:
        if file.endswith('.vcf'):
            valid_vcf(file)
        elif file.endswith('.bed'):
            valid_bed(file)
        else:
            raise Exception('Error: Variant directory has file format(s) other than .vcf and .bed, such as %s' % file)
    label_gtf = 0
    for file in file_dict['label_list']:
        if file.endswith('.gff') or file.endswith('.gff3') or file.endswith('.gtf'):
            label_gtf += 1
            if label_gtf > 1:
                raise Exception('Error: Only a single input gtf/gff file is allow.')
            # valid_gff(file)
        elif file.endswith('.bed'):
            valid_bed(file)
        else:
            raise Exception('Error: Annotation directory has file format(s) other than .gff, .gtf and .bed, such as %s' % file)
    for file in file_dict['filter_list']:
        if file.endswith('.bed'):
            valid_bed(file)
        else:
            raise Exception('Error: Filter directory has file format(s) other than .bed, such as %s' % file)
    return file_dict


def valid_vcf(file):
    with open(file) as f:
        for line in f:
            first = line
            _ = next(f)
            third = next(f)
            vcf_version = float(first.split('v')[1])
            if vcf_version != 4.2:  # fileformat=VCFv4.2
                raise Exception('Error: VCF file %s version is not 4.2' % file)
            nv_version = third.split('-')[1]
            if StrictVersion('1.3.6') > StrictVersion(nv_version):  # source=NanoVar-1.3.6
                raise Exception('Error: VCF file %s was generated by an outdated NanoVar version, version needs to be >= 1.3.6' %
                                file)
            break


def valid_bed(file):
    n = 0
    with open(file) as f:
        for line in f:
            n += 1
            if not line.split('\t')[0].startswith('chr'):
                raise Exception('Error: Contig names has to start with "chr" but detected %s, line %s of %s' % (line.split(
                    '\t')[0], str(n), file))
            try:
                start = int(line.split('\t')[1])
                end = int(line.split('\t')[2])
            except TypeError:
                raise Exception('Error: Non-integer value at line %s of %s' % (str(n), file))
            if start > end:
                raise Exception('Error: End value greater than Start value at line %s of %s' % (str(n), file))


# def valid_gff(file):
#     return None
