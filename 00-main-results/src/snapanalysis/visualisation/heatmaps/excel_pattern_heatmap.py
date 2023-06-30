import numpy as np

from snapanalysis.preprocessing.pulldown_metadata import OUTPUT_FILE as META_FILE
from snapanalysis.visualisation.heatmaps.pattern_heatmap import PULLDOWN_ORDER

import pandas as pd
import itertools
from seaborn.utils import relative_luminance
import seaborn as sns


def create_formats_for_nucleosomes(workbook):
    """
    Creates excel formats for the header that we will use

    :param workbook:
    :return:
    """
    color_palette = pd.read_hdf(META_FILE, '/meta/color_palette')

    # These formats will be used for the color palette
    color_formats = {}
    for key, bg_color in color_palette.items():
        # I love this hack from Seaborn:
        font_color = '#000000' if relative_luminance(bg_color) > .408 else '#FFFFFF'
        color_formats[key] = workbook.add_format(dict(font_color=font_color,
                                                      shrink=True,
                                                      align='center',
                                                      valign='vcenter',
                                                      bg_color=bg_color))
    return color_formats


def _simplify_direction(direction_value):
    direction_value = direction_value.lower()

    if 'forward' in direction_value and 'reverse' not in direction_value:
        return 'Forward'
    elif 'forward' not in direction_value and 'reverse' in direction_value:
        return 'Reverse'
    else:
        raise ValueError('Cannot parse direction from {!r}'.format(direction_value))

def reindex_dataframe_columns(data):
    data = data.copy()
    # Most likely columns will be a multiindex if Direction and Pull-Down ID
    if isinstance(data.columns, pd.MultiIndex):
        if len(data.columns.names) != 2:
            raise ValueError('Can only deal with Direction and Pull-Down ID. Sorry.')

        if 'Direction' not in  data.columns.names or 'Pull-Down ID' not in  data.columns.names:
            # If only two columns, and we have the name for Pull-Down ID, we can infer which column is direciton
            if 'Pull-Down ID' in  data.columns.names:
                if data.columns.names[0] == 'Pull-Down ID':
                    data.columns.names = ['Pull-Down ID', 'Direction']
                else:
                    data.columns.names = ['Direction', 'Pull-Down ID']
                    data.columns = data.columns.swaplevel()
            else:
                raise ValueError('Cannot infer which column is direction, which is Pull-Down ID')

        # Now simplify direction
        levels = data.columns.levels[data.columns.names.index('Direction')]
        data.columns = data.columns.set_levels([_simplify_direction(l) for l in levels], level='Direction')
        data = data.sort_index(axis=1)

        #  This is the only way I can figure out how to reorder predictors by the PD orer
        new_order = []
        for p in PULLDOWN_ORDER:
            for (a, b) in data.columns:
                if a == p:
                    new_order.append((a, b))
        data = data[new_order]
        return data
    # If we have only one index, assume we have just list of PDs
    else:
        return data[[p for p in PULLDOWN_ORDER if p in data.columns]]

def excel_predictors_header(
    reindexed_columns, writer, sheet_name,
    start_row=0, start_col=0, write_direction=True,
    format_ptm_header=None
):
    """
    Draws the header with coloured predictors in excel.

    :param reindexed_columns: columns (order) of the predictors to write, see `reindex_dataframe_columns`
    :param writer: excel writer object
    :param sheet_name: sheet to write to
    :param start_row: row to draw the heatmap
    :param start_col: column to draw the header.
    :param write_direction: If false won't write the header for 'Direction'
    :param format_ptm_header: Formatting for the PTM headers
    :return: number of the row after the header (one to draw data in)
    """
    workbook = writer.book

    # --- Column formats ----
    # For PD ID's etc.
    header_column_format = workbook.add_format({'bold': True, 'bottom': 1,
                                      'align': 'center',
                                      'valign': 'vcenter',
                                      'shrink': True})

    # For nucleosomes
    nucleosome_ptm_formats = create_formats_for_nucleosomes(workbook)

    # --- Read predictors ---

    predictors = pd.read_hdf(META_FILE, '/meta/predictors_categorical_directionalised')
    # Some sanity here
    predictors['Direction'] = predictors['Direction'].replace('Forward', 'F')
    predictors['Direction'] = predictors['Direction'].replace('Reverse', 'R')
    predictors['DNA Methylation'] = predictors['DNA Methylation'].replace('yes', 'm5C')
    predictors['H2A.Z'] = predictors['H2A.Z'].replace('yes', 'H2A.Z')
    predictors = predictors.swaplevel().sort_index()

    standard_predictor_columns = predictors.columns

    # -- Reindex the predictors to match the new header

    # Most likely we will have a multiIndex which has a Pull-Down ID and direction
    if isinstance(reindexed_columns, pd.MultiIndex):
        standard_predictor_columns = predictors.columns
        columns_df = reindexed_columns.to_frame()
        del columns_df['Direction']
        del columns_df['Pull-Down ID']

        columns_df = columns_df.join(predictors, on=['Pull-Down ID', 'Direction'])
        predictors = columns_df[standard_predictor_columns]
        predictors.index = predictors.index.reorder_levels(['Pull-Down ID', 'Direction'])

        predictors = predictors.reindex(reindexed_columns)
    # If we have only one index, assume we have just list of PDs
    else:
        predictors = predictors.droplevel('Direction')
        del predictors['Direction']
        predictors = predictors.drop_duplicates()
        predictors = predictors.reindex(reindexed_columns)

    predictors = predictors.T
    assert reindexed_columns.equals(predictors.columns)

    if sheet_name not in writer.sheets:
        worksheet = workbook.add_worksheet(sheet_name)
        writer.sheets[sheet_name] = worksheet
    else:
        worksheet = writer.sheets[sheet_name]

    for row_i, (row_name, row) in enumerate(predictors.iterrows(), start=start_row+1):

        for (value, pull_down), group_ in itertools.groupby(enumerate(zip(row.values,
                                                                          row.index.get_level_values('Pull-Down ID')),
                                                                      start=start_col+1),
                                                            key=lambda x: x[1]):

            group_ = list(group_)

            start_column = group_[0][0]
            end_column = group_[-1][0]

            if value in ['', 'F', 'R']:
                format_ = nucleosome_ptm_formats['background']
            else:
                format_ = nucleosome_ptm_formats[value]

            if start_column == end_column:
                worksheet.write(row_i, start_column, value, format_)
            else:
                worksheet.merge_range(row_i, start_column,
                                      row_i, end_column,
                                      value,
                                      format_)
    # Headers (PD IDs)

    for col, group_ in itertools.groupby(
            enumerate(predictors.columns.get_level_values('Pull-Down ID'), start=start_col+1),
            key=lambda x: x[1]):

        group_ = list(group_)
        start_column = group_[0][0]
        end_column = group_[-1][0]

        if start_column == end_column:
            worksheet.write(start_row, start_column, col, header_column_format)
        else:
            worksheet.merge_range(start_row, start_column,
                                  start_row, end_column,
                                  col,
                                  header_column_format)

    # Headers (PTMs)

    if format_ptm_header is None:
        header_left = workbook.add_format({'bold': True,
                                           'shrink': True, 'valign': 'vcenter',
                                           'align': 'right'})
    else:
        header_left = format_ptm_header

    for i, row_name in enumerate(predictors.index,
                                 start=start_row+1):
        if row_name == 'Direction' and not write_direction:
            continue

        worksheet.write(i, start_col,
                        row_name,
                        header_left)

    return len(predictors) + start_row + 1, predictors

def excel_heatmap(data,
                  writer, sheet_name,
                  vmin=None, vmid=None, vmax=None,
                  cmap='RdBu_r',
                  start_row=0, start_col=0):

    # Change the order of data columns for what we need
    data = reindex_dataframe_columns(data)

    if isinstance(data.index, pd.MultiIndex):
        index_columns = len(data.index.names)
    else:
        index_columns = 1

    workbook = writer.book
    if sheet_name not in writer.sheets:
        worksheet = workbook.add_worksheet(sheet_name)
        writer.sheets[sheet_name] = worksheet
    else:
        worksheet = writer.sheets[sheet_name]

    first_col = start_col + index_columns
    first_row, __ = excel_predictors_header(data.columns,
                                            writer=writer,
                                            sheet_name=sheet_name,
                                            start_col=first_col-1,
                                            start_row=start_row)

    data.to_excel(writer,
                  sheet_name=sheet_name, startrow=first_row -1, startcol=start_col,
                  header=False, index=True)
    last_row = first_row + len(data) - 1
    last_col = first_col + len(data.columns) - 1

    if cmap is not None:
        # Makes darker colours this way...
        palette = sns.color_palette(cmap, 9).as_hex()
        color_min = palette[0]
        color_mid = palette[4]
        color_max = palette[8]

        if vmin is None:
            vmin = np.min(data.values)
        if vmid is None:
            vmid = np.median(data.values)
        if vmax is None:
            vmax = np.max(data.values)

        worksheet.conditional_format(first_row, first_col, last_row, last_col,
                                     {'type': '3_color_scale',
                                      'min_type': 'num',
                                      'mid_type': 'num',
                                      'max_type': 'num',
                                      'min_value': vmin,
                                      'max_value': vmax,
                                      'mid_value': vmid,
                                      'min_color': color_min,
                                      'mid_color': color_mid,
                                      'max_color': color_max})

    return first_row, first_col, last_row, last_col


