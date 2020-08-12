"""
Scattermaps using seaborn.

A lot of this code is mirrorred in https://gist.github.com/lukauskas/f2f43aad6078a8b5d71b986174487b8c
"""

from seaborn.matrix import _HeatMapper, ClusterGrid
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

from seaborn.external.six import string_types
from seaborn.utils import despine, axis_ticklabels_overlap
import seaborn as sns


class _ScatterMapper(_HeatMapper):
    """
    Draw a scattermap plot, similar to heatmap plot, but use scatter dots instead of heatmap
    """

    def __init__(self, data,
                 marker, marker_size,
                 vmin, vmax, cmap, center, robust, cbar, cbar_kws,
                 xticklabels=True, yticklabels=True, mask=None):

        super(_ScatterMapper, self).__init__(
            data, vmin, vmax, cmap, center, robust, cbar=cbar, cbar_kws=cbar_kws,
            xticklabels=xticklabels, yticklabels=yticklabels, mask=mask,
            # Don't support annotation
            annot=False, fmt=None, annot_kws=None,
        )

        self.marker = marker

        if isinstance(marker_size, float) or isinstance(marker_size, int):
            self.marker_size = marker_size
        elif isinstance(marker_size, pd.DataFrame):
            self.marker_size = marker_size.loc[self.data.index, self.data.columns].values
        else:
            self.marker_size = marker_size

    def plot(self, ax, cax, kws):
        """Draw the scattermap on the provided Axes."""
        # Remove all the Axes spines
        despine(ax=ax, left=True, bottom=True)

        # Draw the heatmap
        data = self.plot_data

        range_y = np.arange(data.shape[0], dtype=int) + 0.5
        range_x = np.arange(data.shape[1], dtype=int) + 0.5
        x, y = np.meshgrid(range_x, range_y)

        hmap = ax.scatter(x, y,
                          c=data,
                          marker=self.marker,
                          cmap=self.cmap,
                          vmin=self.vmin, vmax=self.vmax,
                          s=self.marker_size, **kws)

        # Set the axis limits
        ax.set(xlim=(0, self.data.shape[1]), ylim=(0, self.data.shape[0]))

        # Possibly add a colorbar
        if self.cbar:
            cb = ax.figure.colorbar(hmap, cax, ax, **self.cbar_kws)
            cb.outline.set_linewidth(0)
            # If rasterized is passed to pcolormesh, also rasterize the
            # colorbar to avoid white lines on the PDF rendering
            if kws.get('rasterized', False):
                cb.solids.set_rasterized(True)

        # Add row and column labels
        if isinstance(self.xticks, string_types) and self.xticks == "auto":
            xticks, xticklabels = self._auto_ticks(ax, self.xticklabels, 0)
        else:
            xticks, xticklabels = self.xticks, self.xticklabels

        if isinstance(self.yticks, string_types) and self.yticks == "auto":
            yticks, yticklabels = self._auto_ticks(ax, self.yticklabels, 1)
        else:
            yticks, yticklabels = self.yticks, self.yticklabels

        ax.set(xticks=xticks, yticks=yticks)
        xtl = ax.set_xticklabels(xticklabels)
        ytl = ax.set_yticklabels(yticklabels, rotation="vertical")

        # Possibly rotate them if they overlap
        ax.figure.draw(ax.figure.canvas.get_renderer())
        if axis_ticklabels_overlap(xtl):
            plt.setp(xtl, rotation="vertical")
        if axis_ticklabels_overlap(ytl):
            plt.setp(ytl, rotation="horizontal")

        # Add the axis labels
        ax.set(xlabel=self.xlabel, ylabel=self.ylabel)

        # Annotate the cells with the formatted values
        if self.annot:
            self._annotate_heatmap(ax, hmap)

        # Invert the y axis to show the plot in matrix form
        ax.invert_yaxis()


def scattermap(data,
               marker='o',
               marker_size=100,
               vmin=None, vmax=None, cmap=None, center=None, robust=False,
               linewidths=0, linecolor="white",
               cbar=True, cbar_kws=None, cbar_ax=None,
               square=False, xticklabels="auto", yticklabels="auto",
               mask=None, ax=None, **kwargs):
    """Plot rectangular data as a color-encoded matrix.
    This function is similar to `sns.heatmap`, as it is an Axes-level function that will draw the
    heatmap into the currently-active Axes if none is provided to the ``ax`` argument.
    The main difference is that instead of drawing an actual heatmap with filled squares,
    this function will use the `plt.scatter` behind the scenes to draw a scatterplot-heatmap.
    The default is set to plot a grid of circles, however this can be changed via `marker`
    parameter.
    Parameters
    ----------
    data : rectangular dataset
        2D dataset that can be coerced into an ndarray. If a Pandas DataFrame
        is provided, the index/column information will be used to label the
        columns and rows.
    marker: string, optional
        Marker to use: any marker that `pyplot.scatter` supports. Defaults to circle.
    marker_size: int or rectangular dataset
        Either an integer to set the marker size of all data points to,
        or a 2D dataset (like in `data`) that sets individual point sizes.
        Defaults to 100.
    vmin, vmax : floats, optional
        Values to anchor the colormap, otherwise they are inferred from the
        data and other keyword arguments.
    cmap : matplotlib colormap name or object, or list of colors, optional
        The mapping from data values to color space. If not provided, the
        default will depend on whether ``center`` is set.
    center : float, optional
        The value at which to center the colormap when plotting divergant data.
        Using this parameter will change the default ``cmap`` if none is
        specified.
    robust : bool, optional
        If True and ``vmin`` or ``vmax`` are absent, the colormap range is
        computed with robust quantiles instead of the extreme values.
    linewidths : float, optional
        Width of the border lines that will surround the markers
    linecolor : color, optional
        Color of the border lines to the markers
    cbar : boolean, optional
        Whether to draw a colorbar.
    cbar_kws : dict of key, value mappings, optional
        Keyword arguments for `fig.colorbar`.
    cbar_ax : matplotlib Axes, optional
        Axes in which to draw the colorbar, otherwise take space from the
        main Axes.
    square : boolean, optional
        If True, set the Axes aspect to "equal" so each cell will be
        square-shaped.
    xticklabels, yticklabels : "auto", bool, list-like, or int, optional
        If True, plot the column names of the dataframe. If False, don't plot
        the column names. If list-like, plot these alternate labels as the
        xticklabels. If an integer, use the column names but plot only every
        n label. If "auto", try to densely plot non-overlapping labels.
    mask : boolean array or DataFrame, optional
        If passed, data will not be shown in cells where ``mask`` is True.
        Cells with missing values are automatically masked.
    ax : matplotlib Axes, optional
        Axes in which to draw the plot, otherwise use the currently-active
        Axes.
    kwargs : other keyword arguments
        All other keyword arguments are passed to ``ax.pcolormesh``.
    Returns
    -------
    ax : matplotlib Axes
        Axes object with the heatmap.
    See also
    --------
    clustermap : Plot a matrix using hierachical clustering to arrange the
                 rows and columns.
    Examples
    --------
    Plot a scattermap for a numpy array:
    .. plot::
        :context: close-figs
        >>> import numpy as np; np.random.seed(0)
        >>> import seaborn as sns; sns.set()
        >>> uniform_data = np.random.rand(10, 12)
        >>> ax = scattermap(uniform_data)
    Draw on white axes
    .. plot::
        :context: close-figs
        >>> uniform_data = np.random.rand(10, 12)
        >>> with sns.axes_style("white"):
        ...     ax = scattermap(uniform_data)
    Change the limits of the scattermap:
    .. plot::
        :context: close-figs
        >>> ax = scattermap(uniform_data, vmin=0, vmax=1)
    Plot a scattermap for data centered on 0 with a diverging colormap:
    .. plot::
        :context: close-figs
        >>> normal_data = np.random.randn(10, 12)
        >>> ax = scattermap(normal_data, center=0)
    Plot a dataframe with meaningful row and column labels:
    .. plot::
        :context: close-figs
        >>> flights = sns.load_dataset("flights")
        >>> flights = flights.pivot("month", "year", "passengers")
        >>> ax = scattermap(flights)
    Add border lines around each glyph:
    .. plot::
        :context: close-figs
        >>> ax = scattermap(flights, linewidths=1, linecolor='black')
    Use a different colormap:
    .. plot::
        :context: close-figs
        >>> ax = scattermap(flights, cmap="YlGnBu")
    Center the colormap at a specific value:
    .. plot::
        :context: close-figs
        >>> ax = scattermap(flights, center=flights.loc["January", 1955])
    Plot every other column label and don't plot row labels:
    .. plot::
        :context: close-figs
        >>> data = np.random.randn(50, 20)
        >>> ax = scattermap(data, xticklabels=2, yticklabels=False)
    Don't draw a colorbar:
    .. plot::
        :context: close-figs
        >>> ax = scattermap(flights, cbar=False)
    Use different axes for the colorbar:
    .. plot::
        :context: close-figs
        >>> grid_kws = {"height_ratios": (.9, .05), "hspace": .3}
        >>> f, (ax, cbar_ax) = plt.subplots(2, gridspec_kw=grid_kws)
        >>> ax = scattermap(flights, ax=ax,
        ...                  cbar_ax=cbar_ax,
        ...                  cbar_kws={"orientation": "horizontal"})
    Use a mask to plot only part of a matrix
    .. plot::
        :context: close-figs
        >>> corr = np.corrcoef(np.random.randn(10, 200))
        >>> mask = np.zeros_like(corr)
        >>> mask[np.triu_indices_from(mask)] = True
        >>> with sns.axes_style("white"):
        ...     ax = scattermap(corr, mask=mask, vmax=.3, square=True)
     Change glyph, plot stars instead of circles
    .. plot::
        :context: close-figs
        >>> ax = scattermap(corr, vmax=.3, square=True, marker='*')
    Plot multiple markers on the same plot
    >>> corr = np.corrcoef(np.random.randn(10, 200))
    >>> mask = np.zeros_like(corr)
    >>> mask[np.triu_indices_from(mask)] = True
    >>> with sns.axes_style("white"):
    ...     ax = scattermap(corr, mask=mask, vmax=.3, square=True)
    ...     ax = scattermap(corr, mask=mask.T, vmax=.3, square=True, ax=ax, marker='*', cbar=False)
    Specify size for points
    .. plot::
        :context: close-figs
    >>> with sns.axes_style("white"):
    ...     ax = scattermap(corr, vmax=.3, square=True, marker_size=np.abs(corr)*300)
    """
    # Initialize the plotter object
    plotter = _ScatterMapper(data,
                             marker, marker_size,
                             vmin, vmax, cmap, center, robust,
                             cbar, cbar_kws, xticklabels,
                             yticklabels, mask)

    # Add the pcolormesh kwargs here
    kwargs["linewidths"] = linewidths
    kwargs["edgecolor"] = linecolor

    # Draw the plot and return the Axes
    if ax is None:
        ax = plt.gca()
    if square:
        ax.set_aspect("equal")
    plotter.plot(ax, cbar_ax, kwargs)
    return ax


class TwoDatasetScatterClusterGrid(ClusterGrid):

    def __init__(self, data_agg, data_first, data_second,
                 marker_first, marker_second,
                 marker_size_first=100, marker_size_second=100,
                 pivot_kws=None, z_score=None, standard_scale=None,
                 figsize=None, row_colors=None, col_colors=None, mask=None,
                 row_colors_as_scattermap=False,
                 col_colors_as_scattermap=False,
                 row_colors_scattermap_marker_size=25,
                 col_colors_scattermap_marker_size=25,
                 **kwargs,
                 ):

        super(TwoDatasetScatterClusterGrid, self).__init__(data_agg, pivot_kws, z_score,
                                                           standard_scale,
                                                           figsize, row_colors, col_colors, mask,
                                                           **kwargs)

        assert data_agg.index.equals(data_first.index)
        assert data_agg.index.equals(data_second.index)

        assert data_agg.columns.equals(data_first.columns)
        assert data_agg.columns.equals(data_second.columns)

        self.data_first = data_first
        self.data_first_2d = self.format_data(data_first, pivot_kws, z_score,
                                              standard_scale)

        self.data_second = data_second
        self.data_second_2d = self.format_data(data_second, pivot_kws, z_score,
                                               standard_scale)

        self.marker_first = marker_first
        self.marker_second = marker_second

        if isinstance(marker_size_first, int):
            self.marker_size_first_2d = marker_size_first
        else:
            self.marker_size_first_2d = self.format_data(marker_size_first, pivot_kws)

        if isinstance(marker_size_second, int):
            self.marker_size_2d = marker_size_second
        else:
            self.marker_size_2d = self.format_data(marker_size_second, pivot_kws)

        self.row_colors_as_scattermap = row_colors_as_scattermap
        self.col_colors_as_scattermap = col_colors_as_scattermap
        self.row_colors_scattermap_marker_size = row_colors_scattermap_marker_size
        self.col_colors_scattermap_marker_size = col_colors_scattermap_marker_size

    def plot_colors(self, xind, yind, **kws):
        """Plots color labels between the dendrogram and the heatmap
        Parameters
        ----------
        heatmap_kws : dict
            Keyword arguments heatmap
        """
        # Remove any custom colormap and centering
        kws = kws.copy()
        kws.pop('cmap', None)
        kws.pop('center', None)
        kws.pop('annot', None)
        kws.pop('vmin', None)
        kws.pop('vmax', None)
        kws.pop('robust', None)
        kws.pop('xticklabels', None)
        kws.pop('yticklabels', None)


        # Plot the row colors
        if self.row_colors is not None:
            matrix, cmap = self.color_list_to_matrix_and_cmap(
                self.row_colors, yind, axis=0)

            # Get row_color labels
            if self.row_color_labels is not None:
                row_color_labels = self.row_color_labels
            else:
                row_color_labels = False

            if not self.row_colors_as_scattermap:
                sns.heatmap(matrix, cmap=cmap, cbar=False, ax=self.ax_row_colors,
                            xticklabels=row_color_labels, yticklabels=False, **kws)
            else:
                scattermap(matrix, cmap=cmap, cbar=False, ax=self.ax_row_colors,
                           marker_size=self.row_colors_scattermap_marker_size,
                           xticklabels=row_color_labels, yticklabels=False, **kws)
                self.ax_row_colors.grid(False)

            # Adjust rotation of labels
            if row_color_labels is not False:
                plt.setp(self.ax_row_colors.get_xticklabels(), rotation=90)
        else:
            despine(self.ax_row_colors, left=True, bottom=True)

        # Plot the column colors
        if self.col_colors is not None:
            matrix, cmap = self.color_list_to_matrix_and_cmap(
                self.col_colors, xind, axis=1)

            # Get col_color labels
            if self.col_color_labels is not None:
                col_color_labels = self.col_color_labels
            else:
                col_color_labels = False

            if not self.col_colors_as_scattermap:
                sns.heatmap(matrix, cmap=cmap, cbar=False, ax=self.ax_col_colors,
                    xticklabels=False, yticklabels=col_color_labels, **kws)
            else:
                scattermap(matrix, cmap=cmap, cbar=False, ax=self.ax_col_colors,
                           marker_size=self.col_colors_scattermap_marker_size,
                           xticklabels=False, yticklabels=col_color_labels, **kws)
                self.ax_col_colors.grid(False)



            # Adjust rotation of labels, place on right side
            if col_color_labels is not False:
                self.ax_col_colors.yaxis.tick_right()
                plt.setp(self.ax_col_colors.get_yticklabels(), rotation=0)
        else:
            despine(self.ax_col_colors, left=True, bottom=True)

    def plot_matrix(self, colorbar_kws, xind, yind, **kws):
        self.data2d = self.data2d.iloc[yind, xind]
        self.data_first_2d = self.data_first_2d.iloc[yind, xind]
        self.data_second_2d = self.data_second_2d.iloc[yind, xind]
        self.mask = self.mask.iloc[yind, xind]

        if isinstance(self.marker_size_2d, pd.DataFrame):
            self.marker_size_2d = self.marker_size_2d.iloc[yind, xind]

        # Try to reorganize specified tick labels, if provided
        xtl = kws.pop("xticklabels", "auto")
        try:
            xtl = np.asarray(xtl)[xind]
        except (TypeError, IndexError):
            pass
        ytl = kws.pop("yticklabels", "auto")
        try:
            ytl = np.asarray(ytl)[yind]
        except (TypeError, IndexError):
            pass

        scattermap(self.data_first_2d,
                   marker=self.marker_first,
                   marker_size=self.marker_size_2d,
                   ax=self.ax_heatmap, cbar_ax=self.cax,
                   cbar_kws=colorbar_kws, mask=self.mask,
                   xticklabels=xtl, yticklabels=ytl, **kws)

        scattermap(self.data_second_2d,
                   marker=self.marker_second,
                   marker_size=self.marker_size_2d,
                   ax=self.ax_heatmap, cbar_ax=self.cax,
                   cbar_kws=colorbar_kws, mask=self.mask,
                   xticklabels=xtl, yticklabels=ytl, **kws)

        ytl = self.ax_heatmap.get_yticklabels()
        ytl_rot = None if not ytl else ytl[0].get_rotation()
        self.ax_heatmap.yaxis.set_ticks_position('right')
        self.ax_heatmap.yaxis.set_label_position('right')
        if ytl_rot is not None:
            ytl = self.ax_heatmap.get_yticklabels()
            plt.setp(ytl, rotation=ytl_rot)


class ScatterClusterGrid(ClusterGrid):

    def __init__(self, data, marker_size=100,
                 pivot_kws=None, z_score=None, standard_scale=None,
                 figsize=None, row_colors=None, col_colors=None, mask=None, **kwargs):

        super(ScatterClusterGrid, self).__init__(data, pivot_kws, z_score, standard_scale,
                                                 figsize, row_colors, col_colors, mask, **kwargs)

        if isinstance(marker_size, int):
            self.marker_size_2d = marker_size
        else:
            self.marker_size_2d = self.format_data(marker_size, pivot_kws)

    def plot_colors(self, xind, yind, **kws):
        kws = kws.copy()
        kws.pop('marker', '')

        super().plot_colors(xind, yind, **kws)

    def plot_matrix(self, colorbar_kws, xind, yind, **kws):
        self.data2d = self.data2d.iloc[yind, xind]
        self.mask = self.mask.iloc[yind, xind]

        if isinstance(self.marker_size_2d, pd.DataFrame):
            self.marker_size_2d = self.marker_size_2d.iloc[yind, xind]

        # Try to reorganize specified tick labels, if provided
        xtl = kws.pop("xticklabels", "auto")
        try:
            xtl = np.asarray(xtl)[xind]
        except (TypeError, IndexError):
            pass
        ytl = kws.pop("yticklabels", "auto")
        try:
            ytl = np.asarray(ytl)[yind]
        except (TypeError, IndexError):
            pass

        scattermap(self.data2d, marker_size=self.marker_size_2d,
                   ax=self.ax_heatmap, cbar_ax=self.cax,
                   cbar_kws=colorbar_kws, mask=self.mask,
                   xticklabels=xtl, yticklabels=ytl, **kws)

        ytl = self.ax_heatmap.get_yticklabels()
        ytl_rot = None if not ytl else ytl[0].get_rotation()
        self.ax_heatmap.yaxis.set_ticks_position('right')
        self.ax_heatmap.yaxis.set_label_position('right')
        if ytl_rot is not None:
            ytl = self.ax_heatmap.get_yticklabels()
            plt.setp(ytl, rotation=ytl_rot)


def scatterclustermap(data,
                      pivot_kws=None, method='average', metric='euclidean',
                      z_score=None, standard_scale=None, figsize=None, cbar_kws=None,
                      row_cluster=True, col_cluster=True,
                      row_linkage=None, col_linkage=None,
                      row_colors=None, col_colors=None, mask=None,
                      marker_size=100,
                      expected_size_dendrogram=1.0,
                      expected_size_colors=0.25,
                      **kwargs):
    """Plot a matrix dataset as a hierarchically-clustered heatmap.
    Parameters
    ----------
    data: 2D array-like
        Rectangular data for clustering. Cannot contain NAs.
    pivot_kws : dict, optional
        If `data` is a tidy dataframe, can provide keyword arguments for
        pivot to create a rectangular dataframe.
    method : str, optional
        Linkage method to use for calculating clusters.
        See scipy.cluster.hierarchy.linkage documentation for more information:
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
    metric : str, optional
        Distance metric to use for the data. See
        scipy.spatial.distance.pdist documentation for more options
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html
        To use different metrics (or methods) for rows and columns, you may
        construct each linkage matrix yourself and provide them as
        {row,col}_linkage.
    z_score : int or None, optional
        Either 0 (rows) or 1 (columns). Whether or not to calculate z-scores
        for the rows or the columns. Z scores are: z = (x - mean)/std, so
        values in each row (column) will get the mean of the row (column)
        subtracted, then divided by the standard deviation of the row (column).
        This ensures that each row (column) has mean of 0 and variance of 1.
    standard_scale : int or None, optional
        Either 0 (rows) or 1 (columns). Whether or not to standardize that
        dimension, meaning for each row or column, subtract the minimum and
        divide each by its maximum.
    figsize: tuple of two ints, optional
        Size of the figure to create.
    cbar_kws : dict, optional
        Keyword arguments to pass to ``cbar_kws`` in ``heatmap``, e.g. to
        add a label to the colorbar.
    {row,col}_cluster : bool, optional
        If True, cluster the {rows, columns}.
    {row,col}_linkage : numpy.array, optional
        Precomputed linkage matrix for the rows or columns. See
        scipy.cluster.hierarchy.linkage for specific formats.
    {row,col}_colors : list-like or pandas DataFrame/Series, optional
        List of colors to label for either the rows or columns. Useful to
        evaluate whether samples within a group are clustered together. Can
        use nested lists or DataFrame for multiple color levels of labeling.
        If given as a DataFrame or Series, labels for the colors are extracted
        from the DataFrames column names or from the name of the Series.
        DataFrame/Series colors are also matched to the data by their
        index, ensuring colors are drawn in the correct order.
    mask : boolean array or DataFrame, optional
        If passed, data will not be shown in cells where ``mask`` is True.
        Cells with missing values are automatically masked. Only used for
        visualizing, not for calculating.
    kwargs : other keyword arguments
        All other keyword arguments are passed to ``sns.heatmap``
    Returns
    -------
    clustergrid : ClusterGrid
        A ClusterGrid instance.
    Notes
    -----
    The returned object has a ``savefig`` method that should be used if you
    want to save the figure object without clipping the dendrograms.
    To access the reordered row indices, use:
    ``clustergrid.dendrogram_row.reordered_ind``
    Column indices, use:
    ``clustergrid.dendrogram_col.reordered_ind``
    Examples
    --------
    Plot a clustered heatmap:
    .. plot::
        :context: close-figs
        >>> import seaborn as sns; sns.set(color_codes=True)
        >>> iris = sns.load_dataset("iris")
        >>> species = iris.pop("species")
        >>> g = sns.clustermap(iris)
    Use a different similarity metric:
    .. plot::
        :context: close-figs
        >>> g = sns.clustermap(iris, metric="correlation")
    Use a different clustering method:
    .. plot::
        :context: close-figs
        >>> g = sns.clustermap(iris, method="single")
    Use a different colormap and ignore outliers in colormap limits:
    .. plot::
        :context: close-figs
        >>> g = sns.clustermap(iris, cmap="mako", robust=True)
    Change the size of the figure:
    .. plot::
        :context: close-figs
        >>> g = sns.clustermap(iris, figsize=(6, 7))
    Plot one of the axes in its original organization:
    .. plot::
        :context: close-figs
        >>> g = sns.clustermap(iris, col_cluster=False)
    Add colored labels:
    .. plot::
        :context: close-figs
        >>> lut = dict(zip(species.unique(), "rbg"))
        >>> row_colors = species.map(lut)
        >>> g = sns.clustermap(iris, row_colors=row_colors)
    Standardize the data within the columns:
    .. plot::
        :context: close-figs
        >>> g = sns.clustermap(iris, standard_scale=1)
    Normalize the data within the rows:
    .. plot::
        :context: close-figs
        >>> g = sns.clustermap(iris, z_score=0)
    """

    plotter = ScatterClusterGrid(data, marker_size=marker_size,
                                 pivot_kws=pivot_kws, figsize=figsize,
                                 row_colors=row_colors, col_colors=col_colors,
                                 z_score=z_score, standard_scale=standard_scale,
                                 mask=mask,
                                 expected_size_dendrogram=expected_size_dendrogram,
                                 expected_size_colors=expected_size_colors)

    return plotter.plot(metric=metric, method=method,
                        colorbar_kws=cbar_kws,
                        row_cluster=row_cluster, col_cluster=col_cluster,
                        row_linkage=row_linkage, col_linkage=col_linkage,
                        **kwargs)


def scatterclustermap_two_datasets(data_agg, data_first, data_second,
                                   marker_first, marker_second,
                                   pivot_kws=None, method='average', metric='euclidean',
                                   z_score=None, standard_scale=None, figsize=None, cbar_kws=None,
                                   row_cluster=True, col_cluster=True,
                                   row_linkage=None, col_linkage=None,
                                   row_colors=None, col_colors=None, mask=None,
                                   marker_size_first=100, marker_size_second=100,
                                   expected_size_dendrogram=1.0,
                                   expected_size_colors=0.25,
                                   row_colors_as_scattermap=False,
                                   col_colors_as_scattermap=False,
                                   row_colors_scattermap_marker_size=25,
                                   col_colors_scattermap_marker_size=25,
                                   **kwargs):
    """Plot a matrix dataset as a hierarchically-clustered heatmap.
    Parameters
    ----------
    data_agg: 2D array-like
        Rectangular data for clustering. Cannot contain NAs. Aggregate between data_first and data_second
    data_first: 2D array-like
        Rectangular data for plotting (first).
    data_second: 2D array-like
        Rectangular data for clustering (second).
    marker_first: marker to use for first dataset
    marker_second: marker to use for second dataset
    marker_size_first: marker size to use for first dataset
    marker_size_second: marker size to use for second dataset
    pivot_kws : dict, optional
        If `data` is a tidy dataframe, can provide keyword arguments for
        pivot to create a rectangular dataframe.
    method : str, optional
        Linkage method to use for calculating clusters.
        See scipy.cluster.hierarchy.linkage documentation for more information:
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
    metric : str, optional
        Distance metric to use for the data. See
        scipy.spatial.distance.pdist documentation for more options
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html
        To use different metrics (or methods) for rows and columns, you may
        construct each linkage matrix yourself and provide them as
        {row,col}_linkage.
    z_score : int or None, optional
        Either 0 (rows) or 1 (columns). Whether or not to calculate z-scores
        for the rows or the columns. Z scores are: z = (x - mean)/std, so
        values in each row (column) will get the mean of the row (column)
        subtracted, then divided by the standard deviation of the row (column).
        This ensures that each row (column) has mean of 0 and variance of 1.
    standard_scale : int or None, optional
        Either 0 (rows) or 1 (columns). Whether or not to standardize that
        dimension, meaning for each row or column, subtract the minimum and
        divide each by its maximum.
    figsize: tuple of two ints, optional
        Size of the figure to create.
    cbar_kws : dict, optional
        Keyword arguments to pass to ``cbar_kws`` in ``heatmap``, e.g. to
        add a label to the colorbar.
    {row,col}_cluster : bool, optional
        If True, cluster the {rows, columns}.
    {row,col}_linkage : numpy.array, optional
        Precomputed linkage matrix for the rows or columns. See
        scipy.cluster.hierarchy.linkage for specific formats.
    {row,col}_colors : list-like or pandas DataFrame/Series, optional
        List of colors to label for either the rows or columns. Useful to
        evaluate whether samples within a group are clustered together. Can
        use nested lists or DataFrame for multiple color levels of labeling.
        If given as a DataFrame or Series, labels for the colors are extracted
        from the DataFrames column names or from the name of the Series.
        DataFrame/Series colors are also matched to the data by their
        index, ensuring colors are drawn in the correct order.
    mask : boolean array or DataFrame, optional
        If passed, data will not be shown in cells where ``mask`` is True.
        Cells with missing values are automatically masked. Only used for
        visualizing, not for calculating.
    kwargs : other keyword arguments
        All other keyword arguments are passed to ``sns.heatmap``
    Returns
    -------
    clustergrid : ClusterGrid
        A ClusterGrid instance.
    Notes
    -----
    The returned object has a ``savefig`` method that should be used if you
    want to save the figure object without clipping the dendrograms.
    To access the reordered row indices, use:
    ``clustergrid.dendrogram_row.reordered_ind``
    Column indices, use:
    ``clustergrid.dendrogram_col.reordered_ind``
    Examples
    --------
    Plot a clustered heatmap:
    .. plot::
        :context: close-figs
        >>> import seaborn as sns; sns.set(color_codes=True)
        >>> iris = sns.load_dataset("iris")
        >>> species = iris.pop("species")
        >>> g = sns.clustermap(iris)
    Use a different similarity metric:
    .. plot::
        :context: close-figs
        >>> g = sns.clustermap(iris, metric="correlation")
    Use a different clustering method:
    .. plot::
        :context: close-figs
        >>> g = sns.clustermap(iris, method="single")
    Use a different colormap and ignore outliers in colormap limits:
    .. plot::
        :context: close-figs
        >>> g = sns.clustermap(iris, cmap="mako", robust=True)
    Change the size of the figure:
    .. plot::
        :context: close-figs
        >>> g = sns.clustermap(iris, figsize=(6, 7))
    Plot one of the axes in its original organization:
    .. plot::
        :context: close-figs
        >>> g = sns.clustermap(iris, col_cluster=False)
    Add colored labels:
    .. plot::
        :context: close-figs
        >>> lut = dict(zip(species.unique(), "rbg"))
        >>> row_colors = species.map(lut)
        >>> g = sns.clustermap(iris, row_colors=row_colors)
    Standardize the data within the columns:
    .. plot::
        :context: close-figs
        >>> g = sns.clustermap(iris, standard_scale=1)
    Normalize the data within the rows:
    .. plot::
        :context: close-figs
        >>> g = sns.clustermap(iris, z_score=0)
    """

    plotter = TwoDatasetScatterClusterGrid(data_agg,
                                           data_first,
                                           data_second,
                                           marker_first, marker_second,
                                           marker_size_first=marker_size_first,
                                           marker_size_second=marker_size_second,
                                           pivot_kws=pivot_kws, figsize=figsize,
                                           row_colors=row_colors, col_colors=col_colors,
                                           z_score=z_score, standard_scale=standard_scale,
                                           expected_size_dendrogram=expected_size_dendrogram,
                                           expected_size_colors=expected_size_colors,
                                           row_colors_as_scattermap=row_colors_as_scattermap,
                                           col_colors_as_scattermap=col_colors_as_scattermap,
                                           row_colors_scattermap_marker_size=row_colors_scattermap_marker_size,
                                           col_colors_scattermap_marker_size=col_colors_scattermap_marker_size,
                                           mask=mask)

    return plotter.plot(metric=metric, method=method,
                        colorbar_kws=cbar_kws,
                        row_cluster=row_cluster, col_cluster=col_cluster,
                        row_linkage=row_linkage, col_linkage=col_linkage,
                        **kwargs)
