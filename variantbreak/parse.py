"""
Functions for parsing and annotating variant files

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
import logging
from collections import defaultdict, OrderedDict

from pybedtools import BedTool
import pandas as pd
import numpy as np
import fastcluster


def variant_parse(
        file_dict,
        del_annotate_size,
        merge_buffer,
        annotation_dict,
        gene_info_dict,
        label_col,
        filter_dict,
        filter_col,
        no_annotate_cap,
        cluster_sample,
        auto_filter
):
    """Returns a variant dataframe and a dictionary of sample names to file path.
    """
    print('Preparing variant files')
    logging.info('Preparing variant files')

    # Create sample index list and dictionary
    sample_index_list = []
    sample_index_dict = {}
    sample_name_dict = {}
    for file in file_dict['variant_name_list']:
        index_name = 'S' + str(file_dict['variant_name_list'].index(file) + 1)
        sample_index_list.append(index_name)
        sample_index_dict[file] = index_name

    # Parse variant samples
    data_dict = {}
    merged_bed = None
    bed_class = ''
    for file_path in file_dict['variant_list']:
        n = 0
        sample = sample_index_dict[os.path.basename(file_path)]
        sample_name_dict[sample] = file_path
        if file_path.endswith('.vcf'):
            with open(file_path) as f:
                bed = []
                for line in f:
                    n += 1
                    if not line.startswith('#'):
                        contig1 = contig_rename(line.split('\t')[0])
                        start1 = line.split('\t')[1]
                        sv_id = line.split('\t')[2]
                        score = line.split('\t')[5]
                        info_dict = {}
                        for field in line.split('\t')[7].split(';'):
                            info_dict[field.split('=')[0]] = field.split('=')[1]
                        try:
                            sv_type = info_dict['SVTYPE']
                            end1 = info_dict['END']
                            sv_len = info_dict['SVLEN']
                        except KeyError:
                            logging.critical("Error: Invalid INFO column in %s line %s." % (file_path, n))
                            raise Exception("Error: Invalid INFO column in %s line %s." % (file_path, n))
                        if sv_type == 'BND':
                            try:
                                sv_subtype = info_dict['SV2']
                            except KeyError:
                                logging.critical("Error: Invalid INFO column (SV2) in %s line %s." % (file_path, n))
                                raise Exception("Error: Invalid INFO column (SV2) in %s line %s." % (file_path, n))
                            bnd_record = line.split('\t')[4]
                        else:
                            sv_subtype = '-'
                            bnd_record = '-'
                        try:
                            te = info_dict['TE']
                        except KeyError:
                            te = '-'
                        info_dict['GT'] = line.split('\t')[9].split(':')[0]
                        info_dict['DP'] = line.split('\t')[9].split(':')[1]
                        info_dict['AD'] = line.split('\t')[9].split(':')[2].strip()
                        if sv_type == 'DEL':
                            if int(info_dict['SVLEN'].strip('-')) <= del_annotate_size:  # if del size <= threshold
                                bed.append(contig1 + '\t' + start1 + '\t' + end1 + '\t' + sample + '~' + str(n))
                            else:  # if del size > threshold
                                bed.append(contig1 + '\t' + start1 + '\t' + str(int(start1)+1) + '\t' + sample + '~' + str(n))
                                bed.append(contig1 + '\t' + end1 + '\t' + str(int(end1)+1) + '\t' + sample + '~' + str(n))
                        elif sv_type in ['BND', 'INS']:
                            bed.append(contig1 + '\t' + start1 + '\t' + end1 + '\t' + sample + '~' + str(n))
                        elif sv_type in ['DUP', 'INV']:
                            bed.append(contig1 + '\t' + start1 + '\t' + str(int(start1) + 1) + '\t' + sample + '~' + str(n))
                            bed.append(contig1 + '\t' + end1 + '\t' + str(int(end1) + 1) + '\t' + sample + '~' + str(n))
                        data_dict[sample + '~' + str(n)] = [sample, sv_type, sv_subtype, sv_len, bnd_record,
                                                            score, sv_id, info_dict['GT'], info_dict['DP'],
                                                            info_dict['AD'], te]
            # Create BedTool object
            bed_class = BedTool('\n'.join(bed) + '\n', from_string=True)
        elif file_path.endswith('.bed'):
            with open(file_path) as f:
                bed = []
                n = 0
                for line in f:
                    n += 1
                    if not line.startswith('#'):
                        bed.append('\t'.join(line.split('\t')[0:3]) + '\t' + sample + '~' + str(n))
                    data_dict[sample + '~' + str(n)] = [sample, 'UKN']
            bed_class = BedTool('\n'.join(bed) + '\n', from_string=True)
        if not merged_bed:
            merged_bed = bed_class
        else:  # Merge other BEDs to first BED
            merged_bed = merged_bed.cat(bed_class, d=merge_buffer, c=[4], o="collapse")
    merged_bed = merged_bed.sort()

    # Annotate merged bed
    print('Annotating merged variants')
    logging.info('Annotating merged variants')
    header_dict = {}
    for label in annotation_dict:
        if label == 'GTF':
            for col in ['Gene_id', 'Transcript_id', 'Gene_name', 'Gene_type', 'Gene_feature', 'Feature_no']:
                header_dict[col] = defaultdict(list)
            _intersect = merged_bed.intersect(annotation_dict[label], wa=True, wb=True)
            for line in _intersect:
                line = list(line)
                gene_id = line[7].split('/')[0]
                header_dict['Gene_id'][line[3]].append(gene_id)
                header_dict['Transcript_id'][line[3]].append(line[7].split('/')[1])
                header_dict['Gene_name'][line[3]].append(gene_info_dict[gene_id].split('/')[0] + ';')
                header_dict['Gene_type'][line[3]].append(gene_info_dict[gene_id].split('/')[1])
                header_dict['Gene_feature'][line[3]].append(line[7].split('/')[2])
                header_dict['Feature_no'][line[3]].append(line[7].split('/')[3])
        else:
            header_dict[label] = defaultdict(list)
            _intersect = merged_bed.intersect(annotation_dict[label], wa=True, wb=True)
            for line in _intersect:
                line = list(line)
                header_dict[label][line[3]].append(line[7])

    # Filter merged bed
    print('Filtering merged variants')
    logging.info('Filtering merged variants')
    for filt in filter_dict:
        header_dict[filt] = defaultdict(list)
        _intersect = merged_bed.intersect(filter_dict[filt], wa=True)
        for line in _intersect:
            line = list(line)
            header_dict[filt][line[3]] = ['1']

    # Organize all data into master dataframe
    print('Generating dataframe')
    logging.info('Generating dataframe')
    c = 0
    dummy_list = [[0.0 for _ in range(len(sample_index_list))]]
    df = pd.DataFrame(columns=sample_index_list)
    sv_type_dict = {'DEL': 0.214, 'INV': 0.357, 'INS': 0.500, 'BND': 0.643, 'DUP': 0.786, 'UKN': 0.929}
    header_list = {'Chr': [], 'Start': [], 'End': []}
    for header in label_col + filter_col:
        header_list[header] = []
    for sample in sample_index_list:
        header_list['Hover_' + sample] = []
    for sv in merged_bed:
        c += 1
        sv_name = 'SV' + str(c)
        sv_type = []
        tmp = {}
        _df = pd.DataFrame(dummy_list, index=[sv_name], columns=sample_index_list)
        sv = list(sv)
        for _id in sv[3].split(','):  # If two or more breakends come from the same sample, the last breakend will be included.
            _df.at[sv_name, data_dict[_id.strip()][0]] = sv_type_dict[data_dict[_id.strip()][1]]  # Update internal df
            sv_type.append(data_dict[_id.strip()][1])
            tmp[data_dict[_id.strip()][0]] = data_dict[_id.strip()]
        header_list['Chr'].append(sv[0])
        header_list['Start'].append(sv[1])
        header_list['End'].append(sv[2])
        for header in label_col + filter_col:
            header_list[header].append('/'.join(list(set(header_dict[header][sv[3]]))[0:no_annotate_cap]))
        df = pd.concat([df, _df])  # Update main df
        sv_type = '/'.join(sorted(set(sv_type)))
        rank = 0
        for sample in sample_index_list:
            if _df.loc[sv_name, sample] > 0:
                rank += 1
                annote = []
                for label in annotation_dict:
                    if label == 'GTF':
                        annote.append('/'.join(list(set(header_dict['Gene_name'][sv[3]]))[0:no_annotate_cap]))
                    else:
                        annote.append('/'.join(list(set(header_dict[label][sv[3]]))[0:no_annotate_cap]))
                annote = [g for g in annote if g != '']
                annote = '/'.join(annote).strip(', ').replace(';', '')
                _filter = []
                for filt in filter_dict:
                    if header_dict[filt][sv[3]] == ['1']:
                        _filter.append(filt + ': HIT')
                    else:
                        _filter.append(filt + ': MISS')
                _filter = '; '.join(_filter)
                try:
                    header_list['Hover_' + sample].append('Index: %s<br>SV type: %s<br>Region: %s<br>Annotation: %s<br>'
                                                          '%s<br><br>Sample: %s<br>SV_type: %s<br>SV_subtype: %s<br>'
                                                          'SV_size: %s<br>BND_rec: %s<br>Score: %s<br>SV_id: %s<br>GT: '
                                                          '%s<br>DP: %s<br>AD: %s<br>TE: %s' %
                                                          (sv_name, sv_type, sv[0] + ':' + sv[1] + '-' + sv[2], annote,
                                                           _filter, sample, tmp[sample][1], tmp[sample][2], tmp[sample][3],
                                                           tmp[sample][4], tmp[sample][5], tmp[sample][6], tmp[sample][7],
                                                           tmp[sample][8], tmp[sample][9], tmp[sample][10]))
                except IndexError:
                    header_list['Hover_' + sample].append('Index: %s<br>SV type: %s<br>Region: %s<br>Annotation: %s<br>'
                                                          '%s<br><br>Sample: %s<br>SV_type: %s' %
                                                          (sv_name, sv_type, sv[0] + ':' + sv[1] + '-' + sv[2], annote, _filter,
                                                           sample, tmp[sample][1]))
            else:
                header_list['Hover_' + sample].append('')
        if c % 10000 == 0:
            print('Processed %i merged SV entries' % c)
            logging.info('Processed %i merged SV entries' % c)

    for header in header_list:
        df[header] = header_list[header]

    # Auto Filter out variants intersecting with filter files
    if auto_filter:
        for filt in filter_dict:
            df = df[df[filt] != '1']

    # Obtain clustered reordered sample list
    if cluster_sample:
        # Create numpy array for selected samples
        arr = df.loc[:, sample_index_list].values
        # Change all SV hits to value of 1
        arr[arr != 0] = 1
        # Change astype to int
        arr = arr.astype(int)
        # Transpose array
        arr = np.transpose(arr)
        # Linkage
        z = fastcluster.linkage(arr, 'ward')
        # Get reordered indices
        n = len(z) + 1
        cache = dict()
        for k in range(len(z)):
            c1, c2 = int(z[k][0]), int(z[k][1])
            c1 = [c1] if c1 < n else cache.pop(c1)
            c2 = [c2] if c2 < n else cache.pop(c2)
            cache[n + k] = c1 + c2
        index = cache[2 * len(z)]
        clustered_sample = []
        for i in index:
            clustered_sample.append(sample_index_list[i])
        # Combined new sample order with the rest of column labels
        total_clustered_sample = clustered_sample + [i for i in df.columns if i not in clustered_sample]
    else:
        clustered_sample = sample_index_list
        total_clustered_sample = df.columns.to_list()

    # Reorder sample columns
    df = df.reindex(columns=total_clustered_sample)

    # # Rough sorting main dataframe
    # df.sort_values(by=[sample], inplace=True, ascending=False)

    # Make new copy of df for sample columns
    df_rank = df[clustered_sample].copy()

    # Change all non-zero values to 1
    df_rank[df_rank != 0] = 1

    # Sum all the rows up
    rank = df_rank.sum(axis=1)

    # Add rank list to main df
    df['Rank'] = rank

    # Add sample column binaries to dataframe
    for n in range(len(df_rank.columns)):
        df[str(n)] = df_rank[df_rank.columns[n]]

    # Sort by sample binaries
    for n in range(len(df_rank.columns)):
        df.sort_values(by=[str(n)], ascending=False, inplace=True)

    # Sort by rank
    df.sort_values(by=['Rank'], ascending=False, inplace=True)

    # Convert to 2D array
    arr = df[[str(x) for x in range(len(df_rank.columns))]].values

    # Group similar rows
    c = 1
    cls_dict = {}
    cluster_grp = []
    for i in arr:
        if str(i) not in cls_dict:
            cls_dict[str(i)] = c
            c += 1
        cluster_grp.append(cls_dict[str(i)])

    # Add group list to main df
    df['Cluster'] = cluster_grp

    # Sort by cluster
    df.sort_values(by=["Cluster"], ascending=True, inplace=True)

    # Delete sample binary columns
    df.drop(columns=[str(x) for x in range(len(df_rank.columns))], inplace=True)

    # Replace all nan with blank strings
    df.replace(np.nan, '', inplace=True)

    return df, sample_name_dict


def annotation_parse(
        file_dict,
        promoter_size
):
    """Returns a dictionary of BedTool objects for each annotation file,
    a gene info dictionary, and a label column list.
    """
    print('Preparing annotation files')
    logging.info('Preparing annotation files')
    gene_no, line_holder = 0, 0
    bed, exon_list = [], []
    gene_info_dict = {}
    transcript_start, transcript_end, bed_class = '', '', ''
    annotation_dict = {}
    bed_count = 0
    annote_name_dict = {}
    total_gene_types = set()
    for file_path in file_dict['label_list']:
        if file_path.endswith('.gff') or file_path.endswith('.gff3'):
            file_type = 'gff'
        elif file_path.endswith('.gtf'):
            file_type = 'gtf'
        elif file_path.endswith('.bed'):
            file_type = 'bed'
        else:
            logging.critical('Error: Invalid file extension %s' % file_path)
            raise Exception('Error: Invalid file extension %s' % file_path)
        if file_type in ['gff', 'gtf']:
            with open(file_path) as f:
                n = 0
                for line in f:
                    n += 1
                    if not line.startswith('#'):
                        feature = line.split('\t')[2]
                        if feature == 'gene':
                            gene_no += 1
                            contig = contig_rename(line.split('\t')[0])
                            start = line.split('\t')[3]
                            end = line.split('\t')[4]
                            strand = line.split('\t')[6]
                            info_dict = {}
                            if file_type == 'gtf':
                                for field in filter(None, line.split('\t')[8].strip().split(';')):
                                    info_dict[field.strip().split(' ')[0]] = field.strip().split(' ')[1].strip('\"\'')
                            elif file_type == 'gff':
                                for field in filter(None, line.split('\t')[8].strip().split(';')):
                                    info_dict[field.strip().split('=')[0]] = field.strip().split('=')[1].strip('\"\'')
                            try:
                                gene_id = info_dict['gene_id']
                            except KeyError:
                                logging.critical('Error: File %s is missing gene_id at line %i' % (file_path, n))
                                raise Exception('Error: File %s is missing gene_id at line %i' % (file_path, n))
                            try:
                                gene_type = info_dict['gene_type']
                            except KeyError:
                                try:
                                    gene_type = info_dict['gene_biotype']
                                except KeyError:
                                    logging.critical('Error: File %s is missing gene_type or gene_biotype at line %i' %
                                                     (file_path, n))
                                    raise Exception('Error: File %s is missing gene_type or gene_biotype at line %i' %
                                                    (file_path, n))
                            total_gene_types.add(gene_type)
                            try:
                                gene_name = info_dict['gene_name']
                                if ';' in gene_name:
                                    logging.critical("Error: Gene name %s contains ';' character, please raise an issue on "
                                                     "GitHub." % gene_name)
                                    raise Exception("Error: Gene name %s contains ';' character, please raise an issue on "
                                                    "GitHub." % gene_name)
                            except KeyError:
                                logging.critical('Error: File %s is missing gene_name at line %i' % (file_path, n))
                                raise Exception('Error: File %s is missing gene_name at line %i' % (file_path, n))
                            if strand == '+':
                                bed.append(contig + '\t' + str(max(int(start) - promoter_size, 0)) + '\t' +
                                           str(max(int(start) - 1, 1)) + '\t' + gene_id + '//promoter/')
                            elif strand == '-':
                                bed.append(contig + '\t' + str(int(end) + 1) + '\t' + str(int(end) + promoter_size) + '\t' +
                                           gene_id + '//promoter/')
                            gene_info_dict[gene_id] = gene_name + '/' + gene_type
                        elif feature == 'transcript':
                            transcript_start = line.split('\t')[3]
                            transcript_end = line.split('\t')[4]
                            if file_type == 'gtf':
                                for field in line.split('\t')[8].strip().split(';'):
                                    _field = field.strip().split(' ')[0]
                                    info_dict[_field] = field.strip().split(' ')[1].strip('\"\'')
                                    if _field == 'transcript_id':
                                        break
                            elif file_type == 'gff':
                                for field in line.split('\t')[8].strip().split(';'):
                                    _field = field.strip().split('=')[0]
                                    info_dict[_field] = field.strip().split('=')[1].strip('\"\'')
                                    if _field == 'transcript_id':
                                        break
                            try:
                                transcript_id = info_dict['transcript_id']
                            except KeyError:
                                logging.critical('Error: File %s is missing transcript_id at line %i' % (file_path, n))
                                raise KeyError('Error: File %s is missing transcript_id at line %i' % (file_path, n))
                        elif feature == 'exon':
                            start = line.split('\t')[3]
                            end = line.split('\t')[4]
                            info_dict = {}
                            if file_type == 'gtf':
                                for field in line.split('\t')[8].strip().split(';'):
                                    _field = field.strip().split(' ')[0]
                                    info_dict[_field] = field.strip().split(' ')[1].strip('\"\'')
                                    if _field == 'exon_number':
                                        break
                            elif file_type == 'gff':
                                for field in line.split('\t')[8].strip().split(';'):
                                    _field = field.strip().split('=')[0]
                                    info_dict[_field] = field.strip().split('=')[1].strip('\"\'')
                                    if _field == 'exon_number':
                                        break
                            try:
                                exon_no = info_dict['exon_number']
                            except KeyError:
                                logging.critical('Error: File %s is missing exon_number at line %i' % (file_path, n))
                                raise KeyError('Error: File %s is missing exon_number at line %i' % (file_path, n))
                            bed.append(contig + '\t' + start + '\t' + end + '\t' + gene_id + '/' + transcript_id + '/' + feature
                                       + '/exon' + exon_no)
                            exon_list.append([int(start), int(end)])
                            if strand == '+':
                                if end == transcript_end:
                                    # Generate introns
                                    intron_no = 1
                                    for ind in range(len(exon_list) - 1):
                                        bed.append(contig + '\t' + str(exon_list[ind][1] + 1) + '\t' +
                                                   str(exon_list[ind + 1][0] - 1) + '\t' + gene_id + '/' + transcript_id +
                                                   '/intron/intron' + str(intron_no))
                                        intron_no += 1
                                    exon_list = []
                            elif strand == '-':
                                if start == transcript_start:
                                    # Generate introns
                                    intron_no = 1
                                    for ind in range(len(exon_list) - 1):
                                        bed.append(contig + '\t' + str(exon_list[ind + 1][1] + 1) + '\t' +
                                                   str(exon_list[ind][0] - 1) + '\t' + gene_id + '/' + transcript_id +
                                                   '/intron/intron' + str(intron_no))
                                        intron_no += 1
                                    exon_list = []
                        else:
                            pass
                    if n % 500000 == 0:
                        line_holder += n
                        print('Processed %i GFF entries' % line_holder)
                        logging.info('Processed %i GFF entries' % line_holder)
                        n = 0
            # Create BedTool object
            bed_class = BedTool('\n'.join(bed) + '\n', from_string=True)
            bed_class = bed_class.sort()
            annotation_name = 'GTF'
            # Look for gene types used in VariantMap
            vm_gene_types = {'protein_coding', 'lncRNA', 'miRNA', 'snRNA', 'snoRNA'}
            missing_types = vm_gene_types.difference(total_gene_types)
            if len(missing_types) > 0:
                logging.warning("WARNING: GTF/GFF file missing common gene types such as %s." % ', '.join(list(missing_types)))
        elif file_type == 'bed':
            bed_count += 1
            with open(file_path) as f:
                n = 0
                while n == 0:
                    for line in f:
                        n = 1
                        if line.startswith('#'):
                            try:
                                annotation_name = line.split('\t')[3]  # Use 4th column name to name annotation
                            except IndexError:
                                annotation_name = 'BED' + str(bed_count)
                            if annotation_name.startswith(('GTF', 'BED', 'Filter')):  # If column name is invalid, use default
                                annotation_name = 'BED' + str(bed_count)
                            break
                        else:
                            annotation_name = 'BED' + str(bed_count)
                            break
            with open(file_path) as f:
                n = 0
                for line in f:
                    n += 1
                    if not line.startswith('#'):
                        if len(line.split('\t')) < 4:
                            logging.critical('Error: BED annotation file %s is missing label column on line %s' % (file_path,
                                                                                                                   str(n)))
                            raise Exception('Error: BED annotation file %s is missing label column on line %s' % (file_path,
                                                                                                                  str(n)))
            bed_class = BedTool(file_path)
            bed_class = bed_class.sort()
        annotation_dict[annotation_name] = bed_class
        annote_name_dict[annotation_name] = file_path
    label_col = []
    for name in annotation_dict:
        if name == 'GTF':
            label_col.extend(['Gene_id', 'Transcript_id', 'Gene_name', 'Gene_type', 'Gene_feature', 'Feature_no'])
        else:
            label_col.append(name)
    return annotation_dict, gene_info_dict, label_col, annote_name_dict


def filter_parse(
        file_dict
):
    """Returns a dictionary of BedTool objects for each filter file,
     a list of filter column names, and a dictionary of names to file path.
    """
    print('Preparing filter files')
    logging.info('Preparing filter files')
    filter_count = 0
    filter_dict = OrderedDict()
    filter_col = []
    filter_name_dict = {}
    for file_path in file_dict['filter_list']:
        filter_count += 1
        filter_name = 'Filter' + str(filter_count)
        bed_class = BedTool(file_path)
        bed_class = bed_class.sort()
        filter_dict[filter_name] = bed_class
        filter_name_dict[filter_name] = file_path
    for name in filter_dict:
        filter_col.append(name)
    return filter_dict, filter_col, filter_name_dict


# Attach 'chr' to names of contigs if not present
def contig_rename(name):
    if not str(name).startswith("chr"):
        if int_check(name) or name[0] in ["X", "Y", "M"]:
            return "chr" + str(name)
        else:
            # print('WARNING: Unconventional contig name: %s in %s' % (name, file))
            return name
    else:
        return name


# Check string convertible to int
def int_check(i):
    try:
        int(i)
        return True
    except ValueError:
        return False
