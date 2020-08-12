import pandas as pd
import numpy as np

from snapanalysis.models.ptm_response.main import FDR_THRESHOLD_RESPONSE, FC_THRESHOLD_RESPONSE

from contextlib import contextmanager
from rpy2.robjects.lib import grdevices
from IPython.display import Image, display

from rpy2.robjects import pandas2ri, numpy2ri
from rpy2 import robjects
from rpy2.robjects.packages import importr
import rpy2.robjects.lib.ggplot2 as r_ggplot2
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

r_ggrepel = importr('ggrepel')
r_print = robjects.r['print']
r_grid = importr('grid')
r_base = importr('base')
r_dollar = r_base.__dict__["$"]

r_tick_formatter_code = """
formatterFunOneDigit <- function(x) sprintf("%.1f", x)
formatterFunTwoDigits <- function(x) sprintf("%.2f", x)
"""

r_custom = SignatureTranslatedAnonymousPackage(r_tick_formatter_code, "custom")

@contextmanager
def r_inline_plot(width=4, height=4, dpi=100):

    with grdevices.render_to_bytesio(grdevices.png, 
                                     width=width,
                                     height=height, 
                                     units='in',
                                     res=dpi) as b:

        yield

    data = b.getvalue()
    display(Image(data=data, format='png', embed=True))
    

@contextmanager
def r_plot_to_pdf(filename, width=4, height=4):
    grdevices.pdf(filename, width=width, height=height)
    yield
    grdevices.dev_off()

def transform_data_for_ggplot(data, 
                              label_only_large_fc=False,
                              max_labels=None,
                              special_labels=None,
                              special_labels_mode='all',
                              skip_labels=None,
                              nudges=None):
    data = data.copy()
    
    def _significance_category(row):
        if not row['significant']:
            return 'Background'
        elif not row['significant_and_large_fc']:
            return 'Significant'
        else:
            return 'SignificantFC'
    
    data['significance_category'] = data.apply(_significance_category, axis=1)
    data['abs_logFC'] = data['logFC'].abs()
    data['neg_p_adjust'] = -data['adj.P.Val']
    
    
    data['normed_neg_log10_p_adjust'] = data['neg_log10_p_adjust'] / data['neg_log10_p_adjust'].max()
    data['normed_abs_logFC'] = data['abs_logFC'] / data['abs_logFC'].max()
    
    # Label data by distance from 0,0 coordinate
    data['sort_norm'] = np.sqrt(data['normed_neg_log10_p_adjust'] ** 2 + data['normed_abs_logFC'] ** 2)
    
    
    sorted_data = data.sort_values(by=['significant_and_large_fc', 'significant', 'sort_norm'],
                                      ascending=[False, False, False]).query('significant == True')
    
    if label_only_large_fc: 
        sorted_data = sorted_data.query('significant_and_large_fc == True')    
    
    
    if skip_labels is not None:
        skip_labels = set(skip_labels)
    else:
        skip_labels = set([])
        
    sorted_index = sorted_data.index
    sorted_index = [ix for ix in sorted_index if ix not in skip_labels]
    
    if special_labels_mode != 'as_others' and special_labels is not None:
        sorted_index = [ix for ix in sorted_index if ix not in special_labels]
    
    if max_labels is not None:    
        to_label = sorted_index[:max_labels]
    else:
        to_label = sorted_index
    
    
    to_label = set(to_label)
    if special_labels is not None:
        for label in special_labels:
            if label not in data.index:
                print(f'{label} not in data index..')
                continue
            
            row = data.loc[label]
            
            if special_labels_mode == 'all':
                to_label.add(label)
            elif special_labels_mode == 'significant':
                if row['significant']:
                    to_label.add(label)
            elif special_labels_mode == 'as_others':
                # do not treat special proteins differently
                pass
#         to_label.update(special_labels)
    
        to_label = to_label - skip_labels
    
    data['label'] = ''
    data.loc[to_label, 'label'] = [ix.partition('/')[0] for ix in to_label]
     
    data['group'] = data['significance_category'].copy()
    
    data['nudgex'] = 0
    data['nudgey'] = 0
    
    if nudges is not None:
        for protein, (x,y) in nudges.items():
            data.loc[protein, 'nudgex'] = x
            data.loc[protein, 'nudgey'] = y
    
    if special_labels is not None:
        for protein, group in special_labels.items():
            if protein not in data.index:
                continue
                
            data.loc[protein, 'group'] = group
    
    
    other_groups = ['Background', 'Significant', 'SignificantFC']
    if special_labels is not None:
        special_groups = list(set(special_labels.values()))
    else:
        special_groups = []
        
    group_order = other_groups + special_groups
    
#     print(group_order)
    

    data['sortkey'] = data['group'].apply(group_order.index)
    data = data.sort_values(by=['sortkey', 'abs_logFC'], ascending=[True, False])
    
    with robjects.conversion.localconverter(robjects.default_converter + pandas2ri.converter):
        r_data = pandas2ri.py2rpy(data)

    
    return r_data, data

def r_label_palette(data, palette, color_background='#666666', 
                    color_significant='#B7555A', 
                    color_significant_muted=None):
    
    if palette is None:
        palette = {}
        
    groups = data['group'].unique()
    if color_significant_muted is None:
        color_significant_muted = to_hex(sns.desaturate(color_significant, 0.6))
    palette = palette.copy()
    palette['Background'] = color_background
    palette['Significant'] = color_significant_muted
    palette['SignificantFC'] = color_significant
    
    for group in groups:
        if group not in palette:
            raise Exception(f'No colour given for {group!r}')
            
    return robjects.r.c(**palette)


# --- Volcano plots ----



def plot_volcano_with_r(data, 
                        xlabel='Estimated effect (change in H/L ratio)',
                        title='',
                        max_labels=20,
                        color_background='#737373',
                        color_significant='#252525', 
                        color_significant_muted='#252525',
                        label_only_large_fc=False,
                        special_labels=None,
                        special_palette=None,
                        base_size=12,
                        label_size=3,
                        x='logFC',
                        y='neg_log10_p_adjust',
                        special_labels_mode='all',
                        xlim=None,
                        skip_labels=None,
                        nudges=None,
                        ):
    
 

    
    r_data, r_like_data = transform_data_for_ggplot(data, 
                                                    label_only_large_fc=label_only_large_fc,
                                                    special_labels=special_labels,
                                                    max_labels=max_labels,
                                                    special_labels_mode=special_labels_mode,
                                                    skip_labels=skip_labels,
                                                    nudges=nudges)
    
    plot = r_ggplot2.ggplot(r_data)
    plot += r_ggplot2.theme_minimal(base_size = base_size)
    plot += r_ggplot2.theme(**{'panel.grid.major': r_ggplot2.element_blank(), 
                               'panel.grid.minor': r_ggplot2.element_blank(),
                               'panel.border':  r_ggplot2.element_rect(fill=robjects.rinterface.NA,
                                                                       color = "black")})
    plot += r_ggplot2.theme(text=r_ggplot2.element_text(family='Helvetica', face='plain'))
    plot += r_ggplot2.theme(**{'plot.title': r_ggplot2.element_text(hjust=0.5),
#                               'axis.title.y': r_ggplot2.element_text((t = 0, r = 20, b = 0, l = 0)),
                              })
    
    aes_points = r_ggplot2.aes_string(x=x, y=y, color='group')
    scale_points = r_ggplot2.scale_colour_manual(aes_points, 
                                                 values=r_label_palette(r_like_data, special_palette, 
                                                                        color_background=color_background,
                                                                        color_significant=color_significant,
                                                                        color_significant_muted=color_significant_muted))
    
    plot += aes_points
    plot += scale_points
    
    if xlim is not None:
        plot += r_ggplot2.scale_x_continuous(labels=r_custom.formatterFunTwoDigits, limits=robjects.r.c(*xlim))
    else:
        plot += r_ggplot2.scale_x_continuous(labels=r_custom.formatterFunTwoDigits)
        
    plot += r_ggplot2.scale_y_continuous(labels=r_custom.formatterFunOneDigit)
    
    
    plot += r_ggplot2.geom_hline(yintercept=float(-np.log10(FDR_THRESHOLD_RESPONSE)), 
                                 color='#BDBDBD', alpha=.3)
    plot += r_ggplot2.geom_vline(xintercept=float(FC_THRESHOLD_RESPONSE), 
                             color='#BDBDBD', alpha=.3)
    plot += r_ggplot2.geom_vline(xintercept=-float(FC_THRESHOLD_RESPONSE), 
                             color='#BDBDBD', alpha=.3)
    
    
    plot += r_ggplot2.geom_point(**{'show.legend': False})

    aes_text = r_ggplot2.aes_string(label='label')
    plot += aes_text
    plot += r_ggrepel.geom_text_repel(aes_text, 
                                      nudge_x=r_dollar(r_data, 'nudgex'),
                                      nudge_y=r_dollar(r_data, 'nudgey'),
                                       size=label_size,
                                       family='Helvetica',
                                       **{'show.legend': False, 
                                          'point.padding': 0.25, 
                                          'min.segment.length': 0,
                                          #'max.iter':0,
                                          'segment.color': '#BDBDBD'},
                                     )
    
    plot += r_ggplot2.labs(x=xlabel, y='Adjusted p value (-log10)', title=title)
    
  
        
    plot.plot()
    
def plot_volcano_with_r_inline_and_pdf(filename, *args, width=4, height=4, dpi=100, **kwargs):
    with r_inline_plot(width=width, height=height, dpi=dpi):
        plot_volcano_with_r(*args, **kwargs)
    
    with r_plot_to_pdf(filename, width=width, height=height):
        plot_volcano_with_r(*args, **kwargs)
