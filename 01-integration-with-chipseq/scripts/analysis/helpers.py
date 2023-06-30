import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from tqdm import tqdm
from sklearn.ensemble import IsolationForest
import statsmodels.api as sm
from adjustText import adjust_text
import matplotlib.axes
from scipy.stats import hmean
import pingouin as pg

FIVE_MM_IN_INCH = 0.19685
DPI = 600

MARCS_FEATURE_ORDER = [
    # Same order as in Fig 3
    'H2A.Z', 'meDNA', 
    'H3K4me1', 'H3K4me3', 'H3ac', 'H3K9acK14ac', 'H3K27ac', 
    'H3K9me2', 'H3K9me3', 'H3K27me2', 'H3K27me3',
    'H4ac', 'H4K16ac', 'H4K20me2', 'H4K20me3',         
]


def rotate_by_angle(point, theta):
    # https://scipython.com/book/chapter-6-numpy/examples/creating-a-rotation-matrix-in-numpy/
    c, s = np.cos(theta), np.sin(theta)
    R = np.array(((c, -s), (s, c)))
    
    rotp = R.dot(point)
    return rotp


def make_plot(*, df, 
              x_col, xlabel, y_col, ylabel, 
              hue_col=None, shape_col=None, hue_palette=None, shape_map=None, 
              approx_number_of_labels=10, 
              label_strategy='outliers',
              annot_col=None, title=None, filename=None,
              xticks=None, xticklabels=None,
              lines=None,
              legend=True,
              axes=None,
              do_not_change_limits=False,
              adjust_text_kws=None,
              ):
    
    df = df.dropna(subset=[x_col, y_col])
    
#     outlier_detector = OneClassSVM(nu=outlier_fraction, kernel='rbf', gamma=0.01)
    if approx_number_of_labels >= len(df):
        outliers = set(df.index)
    else:
        if label_strategy == 'outliers':
            outlier_fraction = np.clip(approx_number_of_labels / len(df), 0.01, 0.49)
            outlier_detector = IsolationForest(contamination=outlier_fraction, random_state=42)
            outliers = outlier_detector.fit_predict(df[[x_col, y_col]].values)
            outliers = pd.Series(outliers == -1, index=df.index)
            outliers = set(outliers[outliers].index)
        elif isinstance(label_strategy, list) or isinstance(label_strategy, set):
            outlier_candidates = set(label_strategy) & set(df.index)
            if not outlier_candidates:
                outliers = set()
            else:
                # Get the minimum rank of a label after sorting ascending or descending
                # x_col tends to be more important in most plots, use that value only
                ranks_asc = df.loc[outlier_candidates, [x_col]].rank(ascending=True, pct=True).min(axis=1)
                ranks_desc = df.loc[outlier_candidates, [x_col]].rank(ascending=False, pct=True).min(axis=1)

                min_rank = pd.DataFrame({
                    'asc': ranks_asc,
                    'desc': ranks_desc,
                }).min(axis=1).sort_values()

                outliers = set(min_rank.index[:approx_number_of_labels])
        else:
            raise ValueError(f'{outliers=}')

    # This is useful for the interactive backend of mpl.
    def _fmt_coord(x, y):
        results = df[(df[x_col].between(x-0.05, x+0.05)) & (df[y_col].between(y-0.05, y+0.05))]
        if not results.empty:
            return ','.join(results['Factor_Cell_Identifier_op_b'])
        else:
            return f'x={x:.2f}, y={y:.2f}'
    
        
    xmin, xmax = df[x_col].min(), df[x_col].max()
    xrange = xmax-xmin
    
    ymin, ymax = df[y_col].min(), df[y_col].max()
    yrange = ymax-ymin
    
    if axes is None:
        jgrid = sns.JointGrid(height=10*FIVE_MM_IN_INCH)
        ax_joint = jgrid.ax_joint
        ax_marg_x = jgrid.ax_marg_x
        ax_marg_y = jgrid.ax_marg_y
    elif isinstance(axes, matplotlib.axes.Axes):
        ax_joint = axes
        ax_marg_x = ax_marg_y = None
    else:
        ax_joint, ax_marg_x, ax_marg_y = axes
    
    texts = []
    
    if hue_col is not None:
        groupby_ = df.groupby(hue_col)
        all_hues = df[hue_col].unique()
    else:
        groupby_ = [(None, df)]
        hue_palette = {None: sns.color_palette()[0]}
        all_hues = [None]
        
    if shape_col is None:
        # When shape_col is none, treat shape_map as kwargs
        if shape_map is not None:
            shape_map = {None: shape_map.copy()}
        else:
            shape_map = {None: {}}
    
    if lines is None:
        lines = []
        
        
    for hue_val, subdf in groupby_:
        color = hue_palette[hue_val]
        
        subdf = subdf.sort_values(by=[x_col, y_col])
    
        if ax_marg_x is not None:
            sns.histplot(x=subdf[x_col], ax=ax_marg_x, color=color, alpha=.5)
        if ax_marg_y is not None:
            sns.histplot(y=subdf[y_col], ax=ax_marg_y, color=color, alpha=.5)
        
        if 'lowess' in lines:
            lowess_x, lowess_y = sm.nonparametric.lowess(
                subdf[y_col], subdf[x_col], delta=0.01 * xrange).T
            # This kinda assumes subdf is sorted
            ax_joint.plot(lowess_x, lowess_y, color=color, zorder=2, linewidth=0.75)
        
        if shape_col is None:
            shape_groupby = [(None, subdf)]
        else:
            shape_groupby = subdf.groupby(shape_col)
            
        for sh, subsubdf in shape_groupby:
            kws = shape_map[sh]
            n = len(subsubdf)
            ax_joint.scatter(
                subsubdf[x_col],
                subsubdf[y_col],
                color=color,
                label=f"{hue_val}, {shape_col}={sh} (n={n:,})" if shape_col is not None else f'{hue_val} (n={n:,})',
                rasterized=True,
                **kws,
            )
        
            if annot_col is not None:
                
                for outlier in (outliers & set(subsubdf.index)):
                    text_color = color
                    
                    # Hack in order not to have too bright text for "Other" category
                    if text_color == '#bdbdbd':
                        text_color = 'black'
                    
                    texts.append(
                        ax_joint.text(
                            subsubdf.loc[outlier, x_col], subsubdf.loc[outlier, y_col], subsubdf.loc[outlier, annot_col],
                            color=text_color,
                        )
                    )

    if not do_not_change_limits:
        xlim_min = xmin - xrange * 0.1
        xlim_max = xmax + xrange * 0.1

        ax_joint.set_xlim([xlim_min, xlim_max])

        ylim_min = ymin - yrange * 0.1
        ylim_max = ymax + yrange * 0.1

        ax_joint.set_ylim([ylim_min, ylim_max])
    else:
        xlim_min = xmin
        xlim_max = xmax
        
        ylim_min = ymin
        ylim_max = ymax
    
    
    if texts:
        if adjust_text_kws is None:
            adjust_text_kws = {}
        adjust_text(texts, df[x_col].values, df[y_col].values, ax=ax_joint, 
                    arrowprops=dict(arrowstyle='-', linewidth=0.5, color='#969696'),
                    **adjust_text_kws)

    
   
    if xticks is not None:
        if xticklabels is None:
            xticklabels = xticks
        
        assert len(xticklabels) == len(xticks)
        
        _xticks = []
        _xticklabels = []
        
        # Keep only ticks within limit
        for tick, ticklabel in zip(xticks, xticklabels):
            if xlim_min <= tick <= xlim_max:
                _xticks.append(tick)
                _xticklabels.append(ticklabel)
        
        ax_joint.set_xticks(_xticks)
        ax_joint.set_xticklabels(_xticklabels)
    
    grid_line_kws = dict(color='#969696', linestyle=':', zorder=0, linewidth=0.5)
    
    horizontal_grid_lines = [float(l.partition('=')[2]) for l in lines if l.startswith('y=')]
    for hline in horizontal_grid_lines:
        ax_joint.axhline(hline, **grid_line_kws)
        if ax_marg_y is not None:
            ax_marg_y.axhline(hline, **grid_line_kws)
    
    vertical_grid_lines = [float(l.partition('=')[2]) for l in lines if l.startswith('x=')]
    for vline in vertical_grid_lines:
        ax_joint.axvline(vline, **grid_line_kws)
        if ax_marg_x is not None:
            ax_marg_x.axvline(vline, **grid_line_kws)
    
    if 'xy' in lines:
        ax_joint.plot(
            [np.min([xlim_min, ylim_min]), np.max([xlim_max, ylim_max])], 
            [np.min([xlim_min, ylim_min]), np.max([xlim_max, ylim_max])], 
            **grid_line_kws
        )
        
        ## Shading around diagonal  (positive/negative)
        angle_offset = 0.2 # Same as in other figures in main text
        # Diagonal shading
        absmax = np.max(np.abs([xlim_min, xlim_max, ylim_min, ylim_max]))
        line_drawing_limit = absmax * 2
        vector = np.asarray([1,1]) * line_drawing_limit
        shading_dict = dict(color='#969696', alpha=.1, zorder=0)
            
        for sign_symmetry in [-1.0, 1.0]:
            for sign_angle in [-1.0, 1.0]:
                rotated_vector = rotate_by_angle(sign_symmetry * vector, sign_symmetry * sign_angle * angle_offset)  

                ax_joint.fill_between(
                    [0, rotated_vector[0]],
                    [0, rotated_vector[0]],
                    [0, rotated_vector[1]],
                    **shading_dict,
                )
        ax_joint.set_xlim(xlim_min, xlim_max)
        ax_joint.set_ylim(ylim_min, ylim_max)
            
    if legend:
        handles, labels = ax_joint.get_legend_handles_labels()
        if ax_marg_y:
            ax_marg_y.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
        else:
            ax_joint.legend(handles, labels)
    
    ax_joint.set_xlabel(xlabel)
    ax_joint.set_ylabel(ylabel)
    
    if ax_marg_x:
        ax_marg_x.set_xlabel('')
        for tick in ax_marg_x.get_xticklabels():
            tick.set_visible(False)
            
    if ax_marg_y:
        ax_marg_y.set_ylabel('')
        for tick in ax_marg_y.get_yticklabels():
            tick.set_visible(False)
    
    if title is not None:
        if ax_marg_x:
            ax_marg_x.set_title(title)
        else:
            ax_joint.set_title(title)
    
    if filename is not None:
        print(f"Saving to {filename} and {filename}.tsv")
        plt.savefig(filename, bbox_inches='tight', dpi=DPI)
        df.to_csv(filename + '.tsv', sep='\t')
    
    # useful for interactive backend (%matplotlib widget) only
    ax_joint.format_coord = _fmt_coord
    
    return ax_joint, ax_marg_x, ax_marg_y
    
def null_or_concat(col, sep=';'):
    """
    Concatenates unique values of a string-like column if needed,
    propagates None if null
    """
    
    unique_values = col.unique()
    if pd.isnull(unique_values).all():
        return None
    else:
        return sep.join(map(str, unique_values))
    
def assert_null_or_identical(col, sep=';'):
    """
    Concatenates unique values of a string-like column if needed,
    propagates None if null
    """
    
    unique_values = col.unique()
    if pd.isnull(unique_values).all():
        return None
    else:
        assert len(unique_values) == 1
        return unique_values[0]
    
def nan_aware_hmean(col):
    """
    Computes a NaN aware hmean
    """
    
    col = np.asarray(col)
    col = col[~pd.isnull(col)]
    if len(col) == 0:
        return np.nan
    else:
        return hmean(col)
    
        
def reaggregate_by_factor(df, op=None, factor_col=('metadata', 'factor'), normed_mi_header_name='normalised_mi', special_prefix='marcs'):
    """
    Reaggregates data by `Factor` averaging all ChIP-seq experiments in columns.
    Avoids double counting based on MARCS Feature.
    
    If `op` is None, will default to:
       - harmonic mean for columns that are within `normed_mi_header_name` subheading
       - (normal) mean for other numeric columns
       - For columns that start with `special_prefix` assert that values concatenated are identical or null
       - For other, non-numeric columns ensure concatenation of unique values 
    """
    df = df.copy()
    df.index.name = ('metadata', df.index.name)
    df = df.reset_index()
    
    if op is None:
        
        special_cols = {c for c in df.columns if c[0].startswith(special_prefix)}
        
        numeric_cols = set(df.select_dtypes(include=np.number).columns)
        
        numeric_mi_cols = {c for c in numeric_cols if c[0] == normed_mi_header_name}
        other_numeric_cols = numeric_cols - numeric_mi_cols
        other_non_numeric_cols = set(df.select_dtypes(include=object).columns) - special_cols
        
        # Remove groupby col
        other_non_numeric_cols = other_non_numeric_cols - {factor_col}
        
        op = {}
        found_normed_mi = False
        
        for col in special_cols:
            op[col] = assert_null_or_identical
        for col in numeric_mi_cols:
            op[col] = nan_aware_hmean
        
        for col in other_numeric_cols:
            op[col] = 'mean'
        
        for col in other_non_numeric_cols:
            op[col] = null_or_concat
            
        if not numeric_mi_cols:
            raise Exception(f"Something is wrong op is None, but couldn't find normed mi columns with {normed_mi_header_name!r}")
        
        # For debugging
        print("Aggregations used: ")
        _op_reverse = {}
        for _col, _agg in op.items():
            try:
                _op_reverse[_agg].add(_col)
            except KeyError:
                _op_reverse[_agg] = {_col}
        
        for _agg, _cols in _op_reverse.items():
            print(f'- {_agg} ------------')
            print(f'{_cols}')
    
    _sizes = df.groupby(factor_col).size()
    aggregatable_targets = _sizes[_sizes > 1].index
    print("Only the following {:,} factors have >1 row to aggregate: {}".format(len(aggregatable_targets), sorted(aggregatable_targets)))
    
    return df.groupby(factor_col).agg(op)

def get_default_hue_palette_factor_type():
    palette = {
        'protein': '#bdbdbd',
        'feature_histone': '#4BA27C',
        'feature_accessibility': '#E9B83F',
        'state': '#786D9B',
    }
    return palette

def get_default_shape_palette_factor_type():
    _kws = dict(edgecolor='#CDCDCD', alpha=0.8, linewidth=0.35) 
    return dict(zip(['protein', 'feature_histone', 'feature_accessibility', 'state'], [dict(marker='o', s=4, **_kws), dict(marker='^', s=8, **_kws), dict(marker='s', s=8, **_kws), dict(marker='p', s=8, **_kws)]))


def adjust_pvals(stats_df, *, fdr_method, fdr_alpha):
    """
    Corrects the p-values that are not null
    """
    print(f"Adjusting p-values by {fdr_method=}, {fdr_alpha=}")
    # adjust where values are non-null only
    stats_df = stats_df.copy()
    non_null_p = ~stats_df['p-val'].isnull()
    stats_df.loc[non_null_p, 'significant'], \
        stats_df.loc[non_null_p, 'p-val corrected'] = \
            pg.multicomp(
                stats_df.loc[non_null_p, 'p-val'].values, alpha=fdr_alpha, method=fdr_method,
            )
    
    return stats_df

def get_stats(data, *, column, control_groups, test_groups, alternative='two-sided',
              min_n=5, log2=True, fdr_method='fdr_bh', fdr_alpha=0.05):
    """
        :param data: data aggregated by factor
        :param column: a statistic to summarise (first level column index)
        :param control_groups: the groups to use as control, e.g. ['Neither', 'No data']
        :param test_groups: groups to use as test e.g. [Recruited, Excluded] - note that while control group is pooled, test groups are considered separately
        :param alternative: hypothesis test alternative
        :param min_n: minumum number of observations in the control or test groups
        :param log2: compute log2 differences?
        :param fdr_method: FDR correction method
        :param fdr_alpha: Alpha significance level for the FDR method
    """
    
    # Consider proteins only
    data = data[data[('metadata', 'factor_type')] == 'protein']
    
    stats = []
    
    # Loop by the ChIP factor in columns, in this case this will be a histoen factor
    for chip_factor in data[column].columns:

        # Fetch the scores
        scores = data[(column, chip_factor)]
        mask_scores = ~scores.isnull()
        
        # Loop through MARCS features 
        for marcs_feature in data['marcs_feature_significant_category'].columns:
            
            # Get the "control proteins" i.e. proteins that are in the `control_groups` based on MARCS features
            mask_control = mask_scores & data['marcs_feature_significant_category', marcs_feature].isin(control_groups)
            scores_control = scores[mask_control]
            
            # Calculate some summary statistics
            n_control = len(scores_control)
            mean_control = scores_control.mean()
            std_control = scores_control.std(ddof=1)
            
            # Loop through the test groups (e.g. Recruited, Excluded)
            for group_ in test_groups:
                
                # Get the "test" group of scores
                mask_test = mask_scores & (data['marcs_feature_significant_category', marcs_feature] == group_)
                scores_test = scores[mask_test]
               
                # Get some summary statistics
                n_test = len(scores_test)
                mean_test = scores_test.mean()
                std_test = scores_test.std(ddof=1)
                
                # Prepare the data for saving
                d = {
                    'chip_factor': chip_factor,
                    'marcs_feature': marcs_feature,
                    'group': group_,

                    'n_control': n_control,
                    'mean_control': mean_control,
                    'std_control': std_control,

                    'n_test': n_test,
                    'mean_test': mean_test,
                    'std_test': std_test,
                }

                # If we have reached the minimum number of elements in the group
                if min(n_control, n_test) >= min_n:
                    # Then do the MWU test
                    mwu = pg.mwu(scores_test, scores_control, alternative=alternative).loc['MWU']
                    
                    if log2:
                        print("Computing log2 diffs")
                        log2_scores_test = scores_test.apply(np.log2)
                        log2_scores_control = scores_control.apply(np.log2)
                        d['mean_log2_diff'] = log2_scores_test.mean() - log2_scores_control.mean()
                    else:
                        print("Computing diffs in natural scale")
                        d['mean_diff'] = scores_test.mean() - scores_control.mean()
                    
                    # Add the test statistics to stats
                    for col in mwu.index:
                        d[col] = mwu.loc[col]
                # Store the stats
                stats.append(d)
    
    # At this point the hard work here is done, reindex our stats DF
    stats = pd.DataFrame(stats).set_index(['chip_factor', 'marcs_feature', 'group'])  
    if fdr_method is not None:
        stats = adjust_pvals(stats, fdr_method=fdr_method, fdr_alpha=fdr_alpha)
    
    return stats