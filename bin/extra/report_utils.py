#!/usr/bin/env python3
from collections import Counter

from bokeh.models import ColumnDataSource, FactorRange, Whisker, Range1d, HoverTool
#from bokeh.models import VSpan
import bokeh.transform as bkt
import numpy as np
import pandas as pd
from seaborn.categorical import _BarPlotter

from ezcharts.plots import BokehPlot
from ezcharts.plots import _HistogramPlot, util
from seaborn._statistics import Histogram
from scipy.stats import gaussian_kde
from itertools import cycle
from seaborn.relational import _ScatterPlotter
import numbers

from bokeh.plotting import figure

class MakeRectangles(util.JSCode):
    """
    Make JSCode rectangles.

    Use Instance as an argument to series.renderItem.
    Series.encode should list two x dimensions for bar start and end
    and one y dimension for bar height.
    """

    def __init__(self):
        """Instantiate the class with some jscode.

        Each call to this function creates JS code to make a single echarts
        rectangle.

        First get the x-start, x-end and height of the rectangle in raw data
        coords. eg: var data_x_start = api.value(0);

        Convert the raw data coords to canvas coords.
        Note in canvas coords (x=0, y=0) is in the top left.
        eg: var start = api.coord([data_start, 0]);

        start[1] is y=0 in canvas coords. Subtract the height from this to get
        the canvas coords for the top of the rectangle.
        """
        jscode = """function renderItem(params, api) {
            var data_x_start = api.value(0);
            var data_x_end = api.value(1);
            var data_height = api.value(2);

            var start = api.coord([data_x_start, 0]);
            var end = api.coord([data_x_end, 0]);
            var height = api.size([0, data_height])[1];

            var rectShape = echarts.graphic.clipRectByRect(
                {
                    x: start[0],
                    y: start[1] - height,
                    width: end[0] - start[0],
                    height: height
                },
                {
                    x: params.coordSys.x,
                    y: params.coordSys.y,
                    width: params.coordSys.width,
                    height: params.coordSys.height
                }
            );
            return (
                rectShape && {
                    type: 'rect',
                    transition: ['shape'],
                    shape: rectShape,
                    style: api.style()

                }
            );}"""

        super().__init__(jscode)

def seqkit_stats_boxplot_length(
        data=None, *, x=None, y=None, hue=None, order=None, hue_order=None,
        orient=None, color=None, palette=None, saturation=1, fill=True,
        dodge=None, width=None, gap=None, whis=1.5, linecolor='auto', linewidth=1,
        fliersize=6, hue_norm=None, native_scale=None, log_scale=None,
        formatter=None, legend=None, ax=None, x_title=None, y_title=None, plt_title=None, x_rotation=None, **kwargs):
    """Draw a box plot to show distributions with respect to categories."""

    # deal with stuff we haven't implemented, yet
    not_implemented = [
        hue, hue_order, orient, dodge, width, ax, formatter, legend,
        gap, hue_norm, native_scale, log_scale]
    for i in not_implemented:
        if i is not None:
            raise NotImplementedError(
                f"The parameter with the value {i} is not yet implemented.")

    # use our default palette if no colour options were provided
    if palette is None and color is None:
        palette = BokehPlot.colors

    # we are going to group by our x
    groups = data[x].unique()

    # compute quantiles
    df = data
    # if the user specifies an order for their categorical variables
    if isinstance(order, list):
        # check to make sure user has included all categorical variables
        check = Counter(order) == Counter(groups)
        if check:
            groups = order
        else:
            raise ValueError(
                f"all categories ({groups}) not present in order ({order})")

    # compute IQR outlier bounds
    lower_iqr = df.Q2 - df.Q1
    upper_iqr = df.Q3 - df.Q1

    # whiskers
    if isinstance(whis, float):
        df["upper"] = df.Q3 + whis * upper_iqr
        df["lower"] = df.Q1 - whis * lower_iqr
    else:
        raise NotImplementedError(
            f"'whis' must be float, {type(whis)} is not supported or implemented")

    source = ColumnDataSource(df)

    # quantile boxes color
    if not color:
        color = bkt.factor_cmap(x, palette, groups)

    # use seaborn bits and bobs for color
    if not fill:
        color = None
        linecolor = "black"
    else:
        if linecolor == 'auto':
            linecolor = "black"

    # outlier range
    whisker = Whisker(
        base=x, upper="upper", lower="lower", source=source,
        line_width=linewidth, line_color=linecolor)

    whisker.upper_head.size = whisker.lower_head.size = 20

    y_max = df["upper"].max()
    y_min = df["lower"].min()

    y_range = (y_min, 1.1*y_max)
    #y_range = (0, 25000)

    plt = BokehPlot(x_range=groups, y_range=y_range)
    p = plt._fig
    p.add_layout(whisker)

    p.vbar(
        x, 0.7, "Q2", "Q3", source=source, fill_color=color, line_color=linecolor,
        line_width=linewidth, fill_alpha=saturation)

    p.vbar(
        x, 0.7, "Q1", "Q2", source=source, fill_color=color, line_color=linecolor,
        line_width=linewidth, fill_alpha=saturation)

    p.xgrid.grid_line_color = None
    p.xaxis.axis_label = x.capitalize()

    if plt_title is not None:
        p.title.text = plt_title
    if x_title is not None:
        p.xaxis.axis_label = x_title
    if y_title is not None:
        p.yaxis.axis_label = y_title
    if x_rotation is not None:
        p.xaxis.major_label_orientation = x_rotation
    return plt

def quality_boxplot_from_quantiles(
        data=None, *, x=None, y=None, hue=None, order=None, hue_order=None,
        orient=None, color=None, palette=None, saturation=1, fill=True,
        dodge=None, width=None, gap=None, whis=1.5, linecolor='auto', linewidth=1,
        fliersize=6, hue_norm=None, native_scale=None, log_scale=None,
        formatter=None, legend=None, ax=None, x_title=None, y_title=None, plt_title=None, x_rotation=None, **kwargs):
    """Draw a box plot to show distributions with respect to categories."""

    # deal with stuff we haven't implemented, yet
    not_implemented = [
        hue, hue_order, orient, dodge, width, ax, formatter, legend,
        gap, hue_norm, native_scale, log_scale]
    for i in not_implemented:
        if i is not None:
            raise NotImplementedError(
                f"The parameter with the value {i} is not yet implemented.")

    # use our default palette if no colour options were provided
    if palette is None and color is None:
        palette = BokehPlot.colors

    # we are going to group by our x
    groups = data[x].unique()

    # compute quantiles
    df = data
    # if the user specifies an order for their categorical variables
    if isinstance(order, list):
        # check to make sure user has included all categorical variables
        check = Counter(order) == Counter(groups)
        if check:
            groups = order
        else:
            raise ValueError(
                f"all categories ({groups}) not present in order ({order})")

    # compute IQR outlier bounds
    lower_iqr = df.Q2_qual - df.Q1_qual
    upper_iqr = df.Q3_qual - df.Q1_qual

    # whiskers
    if isinstance(whis, float):
        df["upper"] = df.Q3_qual + whis * upper_iqr
        df["lower"] = df.Q1_qual - whis * lower_iqr
    else:
        raise NotImplementedError(
            f"'whis' must be float, {type(whis)} is not supported or implemented")

    source = ColumnDataSource(df)

    # quantile boxes color
    if not color:
        color = bkt.factor_cmap(x, palette, groups)

    # use seaborn bits and bobs for color
    if not fill:
        color = None
        linecolor = "black"
    else:
        if linecolor == 'auto':
            linecolor = "black"

    # outlier range
    whisker = Whisker(
        base=x, upper="upper", lower="lower", source=source,
        line_width=linewidth, line_color=linecolor)

    whisker.upper_head.size = whisker.lower_head.size = 20

    y_max = df["upper"].max()
    y_min = df["lower"].min()

    y_range = (y_min, 1.1*y_max)
    #y_range = (0, 25000)

    plt = BokehPlot(x_range=groups, y_range=y_range)
    p = plt._fig
    p.add_layout(whisker)

    p.vbar(
        x, 0.7, "Q2_qual", "Q3_qual", source=source, fill_color=color, line_color=linecolor,
        line_width=linewidth, fill_alpha=saturation)

    p.vbar(
        x, 0.7, "Q1_qual", "Q2_qual", source=source, fill_color=color, line_color=linecolor,
        line_width=linewidth, fill_alpha=saturation)

    p.xgrid.grid_line_color = None
    p.xaxis.axis_label = x.capitalize()

    if plt_title is not None:
        p.title.text = plt_title
    if x_title is not None:
        p.xaxis.axis_label = x_title
    if y_title is not None:
        p.yaxis.axis_label = y_title
    if x_rotation is not None:
        p.xaxis.major_label_orientation = x_rotation
    return plt

def create_boxplot(df, column_name, sample, plt_title=None, x_title=None, y_title=None):

    """Create a boxplot for the given column."""

    series = df[column_name]

    plt = BokehPlot(tools="save", x_range=[column_name])
    p = plt._fig
    p.y_range = Range1d(start=series.min() - 10, end=series.max() + 10)

    kde = gaussian_kde(series)
    kde_vals = kde.evaluate(np.linspace(series.min(), series.max(), 100))

    kde_scale = max(kde_vals) / 0.4
    kde_x_pos = 0.5

    kde_vals_left = -kde_vals / kde_scale + kde_x_pos
    kde_vals_right = kde_vals[::-1] / kde_scale + kde_x_pos
    kde_vals = np.hstack([kde_vals_left, kde_vals_right])
    kde_support = np.hstack([
        np.linspace(series.min(), series.max(), 100),
        np.linspace(series.max(), series.min(), 100)
    ])

    q1, q2, q3 = series.quantile([0.25, 0.5, 0.75])
    upper_iqr = q3 - q2 * 2
    lower_iqr = q2 - q1 * 2
    qmin, q1, q2, q3, qmax = series.quantile([0, 0.25, 0.5, 0.75, 1])
    upper = min(qmax, q3 + 1.5 * upper_iqr)
    lower = max(qmin, q1 - 1.5 * lower_iqr)

    source = ColumnDataSource(data={'x': kde_vals, 'y': kde_support})
    
    p.patch('x', 'y', source=source, alpha=0.3)

    padding_top = 10
    p.y_range = Range1d(
            start=series.min() - padding_top,
            end=series.max() + padding_top
    )

    hbar_height = (qmax - qmin) / 500
    whisker_width = 0.1

    p.rect([column_name], lower, whisker_width, hbar_height, line_color="grey")
    p.rect([column_name], upper, whisker_width, hbar_height, line_color="grey")
    p.segment([column_name], upper, [column_name], q3, line_color="grey")
    p.segment([column_name], lower, [column_name], q1, line_color="grey")
    p.vbar([column_name], 0.2, q2, q3, line_color="black")
    p.vbar([column_name], 0.2, q1, q2, line_color="black")

    p.xaxis.major_label_orientation = "vertical"
    p.yaxis.axis_label = 'Values'
    p.xaxis.major_label_text_font_size = "0pt"

    if plt_title is not None:
        p.title.text = plt_title
    if x_title is not None:
        p.xaxis.axis_label = x_title
    if y_title is not None:
        p.yaxis.axis_label = y_title

    return plt



def barplot(
    data=None, *, x=None, y=None, hue=None, order=None, hue_order=None,
    estimator='mean', errorbar=('ci', 95), n_boot=1000, units=None, seed=None,
    orient=None, color=None, palette=None, saturation=1.0, width=0.8,
    errcolor='.26', errwidth=None, capsize=None, dodge=True, ci='deprecated',
    ax=None, nested_x=False, plt_title=None, y_title=None, x_title=None, x_rotation=None, **kwargs,
):
    """Show point estimates as rectangular bars.

    Contrary to the seaborn implementation, setting `dodge=False` does not
    result in overlaying the bars, but rather stacking them.

    nested_x: create a nested plot instead of a grouped plot, which has two X
    axis grouping (hue as outer group)
    """
    # use our default palette if no colour options were provided
    if palette is None and color is None:
        palette = BokehPlot.colors

    # A plot can either be nested or stacked, not both
    if nested_x:
        if hue is None or not dodge:
            raise ValueError(
                "`hue` and `dodge` need to be set when passing `nested_x=True`"
            )
        if orient == "h":
            raise ValueError(
                "`nested_x=True` can only work with `orient='v'`"
            )

    # Create bar plot with seaborn
    plotter = _BarPlotter(
        x,
        y,
        hue,
        data,
        order,
        hue_order,
        estimator,
        errorbar,
        n_boot,
        units,
        seed,
        orient,
        color,
        palette,
        saturation,
        width,
        errcolor,
        errwidth,
        capsize,
        dodge,
    )

    # Nested X labeling requires the factors to be provided as a list of
    # tuples with the hue and x levels for each value in the dataframe:
    # factors = [
    #   ("hue1", "x1"),
    #   ("hue1", "x2"),
    #   ("hue2", "x1"),
    #   ("hue2", "x2"),
    # ]
    if nested_x:
        factors = [
            (str(x0), str(x1)) for x0 in plotter.hue_names for x1 in plotter.group_names
        ]
        group_names = FactorRange(*factors)
    else:
        group_names = [str(x) for x in plotter.group_names]

    if plotter.orient == "v":
        plt = BokehPlot(x_range=group_names)
    else:
        plt = BokehPlot(y_range=group_names)

    p = plt._fig

    # Define plot orientation
    # If both hue and dodge=False are provided, make a stacked bar chart
    if plotter.orient == "v":
        plot_bars_func = p.vbar_stack if not dodge and hue else p.vbar
    else:
        plot_bars_func = p.hbar_stack if not dodge and hue else p.hbar

    if plotter.plot_hues is None:
        #data = dict(groups=plotter.group_names)
        data = ColumnDataSource(data)
        data.data["groups"]=plotter.group_names
        #print(data.data)
        #print(data.data[x])
        # simple barplot (i.e. only a single group of bars)
        if plotter.orient == "v":
            plot_bars_func(
                x=x,
                top=y,
                source = data,
                width=0.9,
                color=bkt.factor_cmap(
                   "groups", palette=plotter.colors.as_hex(), factors=plotter.group_names, end=1
                ),
                **kwargs
            )
        else:
            plot_bars_func(
                y=group_names,
                right=plotter.statistic,
                height=0.9,
                color=plotter.colors.as_hex(),
                source=data,
                **kwargs
            )
    elif nested_x:
        # Nested X uses different labelling of the input values
        # x: list of factors described above
        # y: values in a single list/vector
        # legend_name: point to to the legend labelling
        # Create the legend labels
        dual_legend_label = [
            str(x0) for x0 in plotter.hue_names for _ in plotter.group_names
        ]
        # Create the data structure
        data = dict(
            x=factors, y=np.hstack(plotter.statistic.T), legend_name=dual_legend_label
        )
        # Convert to ColumnDataSource
        data = ColumnDataSource(data)
        plot_bars_func(
            x="x",
            top="y",
            source=data,
            line_color="white",
            legend_field="legend_name",
            color=bkt.factor_cmap(
                "x", palette=plotter.colors.as_hex(), factors=plotter.hue_names, end=1
            ),
            **kwargs,
        )
    else:
        data = dict(groups=plotter.group_names)
        data.update(dict(zip(plotter.hue_names, plotter.statistic.T)))
        data = ColumnDataSource(data)
        #print(data)
        if not dodge:
            # stacked bars (define orientation-specific kwargs first)

            if plotter.orient == "v":
                extra_kwargs = dict(x="groups", width=0.95)
            else:
                extra_kwargs = dict(y="groups", height=0.95)
            plot_bars_func(
                plotter.hue_names,
                source=data,
                color=plotter.colors.as_hex(),
                legend_label=plotter.hue_names,
                **extra_kwargs,
                **kwargs,
            )
        else:
            # grouped bars (we need to add a series for each hue level)
            for hue_level, offset, color in zip(
                plotter.hue_names,
                plotter.hue_offsets,
                plotter.colors.as_hex(),
            ):
                # define orientation-specific kwargs
                if plotter.orient == "v":
                    extra_kwargs = dict(
                        x=bkt.dodge("groups", offset, range=p.x_range),
                        top=hue_level,
                        width=plotter.nested_width,
                    )
                else:
                    extra_kwargs = dict(
                        y=bkt.dodge("groups", offset, range=p.y_range),
                        right=hue_level,
                        height=plotter.nested_width,
                    )
                plot_bars_func(
                    source=data,
                    color=color,
                    line_color="white",
                    legend_label=hue_level,
                    **extra_kwargs,
                    **kwargs,
                )

        p.legend.orientation = "horizontal"
        # even though we call `add_layout()` to pull the legend outside of the plot area
        # below, we still need to set `location` to "top" here to make sure the legend
        # is centered properly
        p.legend.location = "top"
        p.add_layout(p.legend[0], "above")

    if plotter.orient == "v":
        p.xgrid.grid_line_color = None
        p.y_range.start = 0
    else:
        p.ygrid.grid_line_color = None
        p.x_range.start = 0

    if not nested_x:
        p.xaxis.axis_label = x.capitalize()
        p.yaxis.axis_label = y.capitalize()
    
    if plt_title is not None:
        p.title.text = plt_title
    if x_title is not None:
        p.xaxis.axis_label = x_title
    if y_title is not None:
        p.yaxis.axis_label = y_title
    if x_rotation is not None:
        p.xaxis.major_label_orientation = x_rotation

    return plt


# this function can be modified
def boxplot(
        data=None, *, x=None, y=None, hue=None, order=None, hue_order=None,
        orient=None, color=None, palette=None, saturation=1, fill=True,
        dodge=None, width=None, gap=None, whis=1.5, linecolor='auto', linewidth=1,
        fliersize=6, hue_norm=None, native_scale=None, log_scale=None,
        formatter=None, legend=None, ax=None, **kwargs):
    """Draw a box plot to show distributions with respect to categories."""
    # deal with stuff we haven't implemented, yet
    not_implemented = [
        hue, hue_order, orient, dodge, width, ax, formatter, legend,
        gap, hue_norm, native_scale, log_scale]
    for i in not_implemented:
        if i is not None:
            raise NotImplementedError(
                f"The parameter with the value {i} is not yet implemented.")

    # use our default palette if no colour options were provided
    if palette is None and color is None:
        palette = BokehPlot.colors

    # we are going to group by our x
    groups = data[x].unique()

    # compute quantiles
    qs = data.groupby(x)[y].quantile([0.25, 0.5, 0.75])
    qs = qs.unstack().reset_index()
    qs.columns = [x, "q1", "q2", "q3"]
    df = pd.merge(data, qs, on=x, how="left")

    # if the user specifies an order for their categorical variables
    if isinstance(order, list):
        # check to make sure user has included all categorical variables
        check = Counter(order) == Counter(groups)
        if check:
            groups = order
        else:
            raise ValueError(
                f"all categories ({groups}) not present in order ({order})")

    # compute IQR outlier bounds
    iqr = df.q3 - df.q1

    # whiskers
    if isinstance(whis, float):
        df["upper"] = df.q3 + whis * iqr
        df["lower"] = df.q1 - whis * iqr
    else:
        raise NotImplementedError(
            f"'whis' must be float, {type(whis)} is not supported or implemented")

    source = ColumnDataSource(df)

    # quantile boxes color
    if not color:
        color = bkt.factor_cmap(x, palette, groups)

    # use seaborn bits and bobs for color
    if not fill:
        color = None
        linecolor = "black"
    else:
        if linecolor == 'auto':
            linecolor = "black"

    # outlier range
    whisker = Whisker(
        base=x, upper="upper", lower="lower", source=source,
        line_width=linewidth, line_color=linecolor)

    whisker.upper_head.size = whisker.lower_head.size = 20

    # outliers
    outliers = df[~df[y].between(df.lower, df.upper)]

    # we need to accpount for whskers and outliers in the y-range
    y_min = df["lower"].min()
    y_max = df["upper"].max()

    # if there are outliers then we need to take those into account
    if not outliers.empty:
        y_min = np.min([outliers[y].min(), df["lower"].min()])
        y_max = np.max([outliers[y].max(), df["upper"].max()])

    y_range = (1.1 * y_min, 1.1 * y_max)

    plt = BokehPlot(x_range=groups, y_range=y_range)
    p = plt._fig
    p.add_layout(whisker)

    p.vbar(
        x, 0.7, "q2", "q3", source=source, fill_color=color, line_color=linecolor,
        line_width=linewidth, fill_alpha=saturation)

    p.vbar(
        x, 0.7, "q1", "q2", source=source, fill_color=color, line_color=linecolor,
        line_width=linewidth, fill_alpha=saturation)

    p.scatter(x, y, source=outliers, size=fliersize, color="black", alpha=0.5)

    p.xgrid.grid_line_color = None
    p.xaxis.axis_label = x.capitalize()
    p.yaxis.axis_label = y.capitalize()

    return plt

def telo_barplot(
    data=None, *, x=None, y=None, hue=None, order=None, hue_order=None,
    estimator='mean', errorbar=('ci', 95), n_boot=1000, units=None, seed=None,
    orient=None, color=None, palette=None, saturation=1.0, width=0.8,
    errcolor='.26', errwidth=None, capsize=None, dodge=True, ci='deprecated',
    ax=None, nested_x=False, plt_title=None, y_title=None, x_title=None, x_rotation=None,
    legend_loc="top", legend_orientation="horizontal", hide_x_tick_labels=False, **kwargs,
):
    """Show point estimates as rectangular bars.

    Contrary to the seaborn implementation, setting `dodge=False` does not
    result in overlaying the bars, but rather stacking them.

    nested_x: create a nested plot instead of a grouped plot, which has two X
    axis grouping (hue as outer group)
    """
    # use our default palette if no colour options were provided
    if palette is None and color is None:
        palette = BokehPlot.colors

    # Create bar plot with seaborn
    plotter = _BarPlotter(
        x,
        y,
        hue,
        data,
        order,
        hue_order,
        estimator,
        errorbar,
        n_boot,
        units,
        seed,
        orient,
        color,
        palette,
        saturation,
        width,
        errcolor,
        errwidth,
        capsize,
        dodge,
    )

    group_names = [str(x) for x in plotter.group_names]

    plt = BokehPlot(x_range=group_names,tooltips=[("Sample", "@groups"),("Bin Start", "$name"),("% of Telomeres", "@$name")], tools="hover")

    p = plt._fig

    # Define plot orientation
    # If both hue and dodge=False are provided, make a stacked bar chart
    plot_bars_func = p.vbar_stack 

    data = dict(groups=plotter.group_names)
    data.update(dict(zip(plotter.hue_names, plotter.statistic.T)))
    data = ColumnDataSource(data)

    plot_bars_func(
        plotter.hue_names,
        x='groups',
        source=data,
        color=plotter.colors.as_hex(),
        legend_label=plotter.hue_names,
        width=0.95,
        **kwargs,
    )
    
    p.legend.orientation = legend_orientation
    p.legend.location = legend_loc
    if legend_loc == "right":
        p.add_layout(p.legend[0], "right")
    else:
        p.add_layout(p.legend[0], "above")

    p.xgrid.grid_line_color = None
    p.y_range.start = 0

    p.xaxis.axis_label = x.capitalize()
    p.yaxis.axis_label = y.capitalize()
    
    if plt_title is not None:
        p.title.text = plt_title
    if x_title is not None:
        p.xaxis.axis_label = x_title
    if y_title is not None:
        p.yaxis.axis_label = y_title
    if x_rotation is not None:
        p.xaxis.major_label_orientation = x_rotation
    if hide_x_tick_labels:
        p.xaxis.major_label_text_font_size = "0pt"
    p.legend.title = "Telomere Length"
    
    return plt


def telo_length_hist(data=None, *, x=None, y=None, hue=None, weights=None,
    stat='count', bins='auto', binwidth=None, binrange=None,
    discrete=None, cumulative=False, common_bins=True,
    common_norm=True, multiple='layer', element='bars',
    fill=True, shrink=1, kde=False, kde_kws=None,
    line_kws=None, thresh=0, pthresh=None, pmax=None,
    cbar=False, cbar_ax=None, cbar_kws=None, palette=None,
    hue_order=None, hue_norm=None, color=None, log_scale=None,
        legend=True, ax=None, plt_title=None, x_title=None, y_title=None, **quad_kwargs):
    """Plot univariate or multivariate histograms."""
    plt = BokehPlot()

    plot = figure()
    
    # print(plot)
    # print(type(plot))
    # glyph = VSpan(x=np.mean(data))
    # plot.add_glyph(glyph)

    # print(plt._fig)
    # print(type(plt._fig))

    estimate_kws = dict(
        stat=stat,
        bins=bins,
        binwidth=binwidth,
        binrange=binrange,
        discrete=discrete,
        cumulative=cumulative,
    )
    
    data = pd.DataFrame(data)

    estimator = Histogram(**estimate_kws)

    if data.ndim > 1 and data.shape[1] > 1:
        # multivariate data
        opacity = 0.5
        if palette is None:
            palette = util.choose_palette()
    else:
        opacity = 1.0
        if color is None:
            palette = util.choose_palette()
        else:
            palette = [color]
    if hue:
        data = data.pivot(columns=hue, values=data.columns[0])
    # this just looks over values if data is 1D
    # for var, color in zip(data, cycle(palette)):
    for col, color in zip(data.columns, cycle(palette)):
        quad_kwargs = {}
        if len(data.columns) > 1:
            quad_kwargs["legend_label"] = col
        variable_data = data[col].dropna()
        heights, edges = estimator(variable_data, weights=weights)

        plt._fig.quad(
            top=heights, bottom=0, left=edges[:-1], right=edges[1:],
            fill_color=color, fill_alpha=opacity, line_color=color, **quad_kwargs
        )
    
    #plt._fig.vspan(x=[np.mean(data)], line_width=[1], color="red")
    plt._fig.y_range.start = 0
    hover = plt._fig.select(dict(type=HoverTool))
    hover.tooltips = [(stat.capitalize(), "@top")]

    if plt_title is not None:
        plt._fig.title.text = plt_title
    if x_title is not None:
        plt._fig.xaxis.axis_label = x_title
    if y_title is not None:
        plt._fig.yaxis.axis_label = y_title

    return plt

DEFAULT_PALETTE = util.choose_palette()
DEFAULT_COLOR = DEFAULT_PALETTE[0]

class Mixin:
    """Mixin class to define how to plot lines and points."""

    def plot(self, plt: BokehPlot, kws):
        """Plot the plot."""
        # Get the potentially normalised data from seaborn
        data = self.plot_data.dropna()
        if data.empty:
            return

        if 'hue' not in data:
            # add a `hue` column (this can't be set to `None` as it would produce an
            # empty plot)
            data['hue'] = "values"

        x_name = self.variables.get("x", None)
        y_name = self.variables.get("y", None)

        symbol_size = kws.get('s', 8)
        if not isinstance(symbol_size, numbers.Number):
            raise NotImplementedError(
                "symbol size (s) currently only accepts a constant value")

        line_width = kws.get('linewidth', 3)
        if not isinstance(line_width, numbers.Number):
            raise NotImplementedError(
                "linewidth currently only accepts a constant value")

        # Available markers/glyphs can be found here:
        # https://docs.bokeh.org/en/latest/docs/reference/models/glyphs.html
        marker = kws.get('marker')
        # TODO: size and style are also grouping "semantic" variables

        y_min = 0
        for hue_name, df in data.groupby('hue'):
            relational_kwargs = {}
            if data['hue'].nunique() > 1:
                relational_kwargs["legend_label"] = hue_name

            y_min = min(y_min, df.y.min())
            color = kws.get("color")
            if self._hue_map.levels:
                # the hue map is not empty (i.e. it was initialised with a palette); use
                # this to override the color (in case one was passed)
                color = self._hue_map(hue_name)

            if self.series_type == 'line':
                plt._fig.line(
                    df.x, df.y, line_color=color,
                    line_width=line_width, **relational_kwargs)

            elif self.series_type == 'scatter':
                plt._fig.scatter(
                    df.x, df.y, color=color, **relational_kwargs)

            if marker is not False:
                if marker is None or marker is True:
                    marker = 'circle'
                # Use scatter to add the markers
                plt._fig.scatter(
                    df.x, df.y, size=symbol_size, marker=marker,
                    color=color, **relational_kwargs)

            plt._fig.xaxis.axis_label = x_name
            plt._fig.yaxis.axis_label = y_name

        # Default the y_range start to zero unless y contains negative values.
        plt._fig.y_range.start = y_min

        # TODO: add legend
        return plt
    
def scatterplot(
        data=None, *,
        x=None, y=None, hue=None, size=None, style=None,
        palette=None, hue_order=None, hue_norm=None,
        sizes=None, size_order=None, size_norm=None,
        markers=True, style_order=None, legend="auto", ax=None, hover_tooltips=None,
        plt_title=None, x_title=None, y_title=None,
        **kwargs):
    """Draw a scatter plot with possibility of several semantic groupings."""
    # see https://github.com/mwaskom/seaborn/blob/949dec3666ab12a366d2fc05ef18d6e90625b5fa/seaborn/relational.py#L726  # noqa

    if palette is None:
        palette = DEFAULT_PALETTE
    
    if not kwargs.get("color"):
        kwargs["color"] = DEFAULT_COLOR
    
    variables = ScatterPlotter.get_semantics(locals())
    p = ScatterPlotter(data=data, variables=variables, legend=legend)
    

    p.map_hue(palette=palette, order=hue_order, norm=hue_norm)
    p.map_size(sizes=sizes, order=size_order, norm=size_norm)
    p.map_style(markers=markers, order=style_order)
    
    plt = BokehPlot(title=kwargs.get('title', ""))
    hover = plt._fig.select(dict(type=HoverTool))
    
    if hover_tooltips is None:
        hover.tooltips = [(x, "@x"), (y, "@y")]
    else:
        hover.tooltips = hover_tooltips
    
    p.plot(plt, kwargs)

    if plt_title is not None:
        plt._fig.title.text = plt_title
    if x_title is not None:
        plt._fig.xaxis.axis_label = x_title
    if y_title is not None:
        plt._fig.yaxis.axis_label = y_title

    return plt


class ScatterPlotter(Mixin, _ScatterPlotter):
    """Making scatter plots."""

    series_type = 'scatter'



def create_multisample_boxplot(df, column_names, min_q, max_q, 
                                plt_title=None, x_title=None, y_title=None,
                                x_rotation=None, x_labels=None):

    """Create a boxplot for the given column."""


    plt = BokehPlot(tools="save", x_range=column_names)
    p = plt._fig

    p.y_range = Range1d(start=min_q, end=max_q)

    for column in column_names:
        series = df[column]

        q1, q3 = series.quantile([0.25, 0.75])
        iqr = q3 - q1
        qmin, q1, q2, q3, qmax = series.quantile([0, 0.25, 0.5, 0.75, 1])
        upper = min(qmax, q3 + 1.5 * iqr)
        lower = max(qmin, q1 - 1.5 * iqr)

        hbar_height = 0.2
        whisker_width = 0.1

        p.rect([column], lower, whisker_width, hbar_height, line_color="grey")
        p.rect([column], upper, whisker_width, hbar_height, line_color="grey")
        p.segment([column], upper, [column], q3, line_color="grey")
        p.segment([column], lower, [column], q1, line_color="grey")
        p.vbar([column], 0.2, q2, q3, line_color="black")
        p.vbar([column], 0.2, q1, q2, line_color="black")

    #p.xaxis.major_label_overrides = x_labels
    #p.xaxis.major_label_orientation = "vertical"
    #p.xaxis.ticker = list(range(len(x_labels)))

    if plt_title is not None:
        p.title.text = plt_title
    if x_title is not None:
        p.xaxis.axis_label = x_title
    if y_title is not None:
        p.yaxis.axis_label = y_title
    if x_rotation is not None:
        p.xaxis.major_label_orientation = x_rotation

    return plt



def get_vals(series, kde_x_pos):

    kde = gaussian_kde(series)
    kde_vals = kde.evaluate(np.linspace(series.min(), series.max(), 100))

    kde_scale = max(kde_vals) / 0.4

    kde_vals_left = -kde_vals / kde_scale + kde_x_pos
    kde_vals_right = kde_vals[::-1] / kde_scale + kde_x_pos
    kde_vals = np.hstack([kde_vals_left, kde_vals_right])
    kde_support = np.hstack([
        np.linspace(series.min(), series.max(), 100),
        np.linspace(series.max(), series.min(), 100)
    ])

    q1, q2, q3 = series.quantile([0.25, 0.5, 0.75])
    upper_iqr = (q3 - q2) * 2
    lower_iqr = (q2 - q1) * 2
    qmin, q1, q2, q3, qmax = series.quantile([0, 0.25, 0.5, 0.75, 1])
    upper = min(qmax, q3 + (1.5 * upper_iqr))
    lower = max(qmin, q1 - (1.5 * lower_iqr))

    hbar_height = (qmax - qmin) / 500

    return {"kde_vals":kde_vals, "kde_support":kde_support, "lower":lower, "upper":upper, "q1":q1, "q2":q2, "q3":q3, "hbar_height":hbar_height}

def create_boxplot_by_strand(df, column_name, plt_title=None, x_title=None, y_title=None):

    """Create a boxplot for the given column."""


    plt = BokehPlot(tools="save", x_range=["G Strand", "C Strand"])
    p = plt._fig
    p.y_range = Range1d(start=df[column_name].min() - 10, end=df[column_name].max() + 10)

    g_strand = get_vals(df[df["strand"]=="G"][column_name], 0.5)
    c_strand = get_vals(df[df["strand"]=="C"][column_name], 1.5)
    
    p.patch(g_strand["kde_vals"], g_strand["kde_support"], alpha=0.3)
    p.patch(c_strand["kde_vals"], c_strand["kde_support"], alpha=0.3)

    padding_top = 10
    p.y_range = Range1d(
            start=df[column_name].min() - padding_top,
            end=df[column_name].max() + padding_top
    )

    whisker_width = 0.1
    
    p.rect(["G Strand"], g_strand["lower"], whisker_width, g_strand["hbar_height"], line_color="grey")
    p.rect(["G Strand"], g_strand["upper"], whisker_width, g_strand["hbar_height"], line_color="grey")
    p.segment(["G Strand"], g_strand["upper"], ["G Strand"], g_strand["q3"], line_color="grey")
    p.segment(["G Strand"], g_strand["lower"], ["G Strand"], g_strand["q1"], line_color="grey")
    p.vbar(["G Strand"], 0.2, g_strand["q2"], g_strand["q3"], line_color="black")
    p.vbar(["G Strand"], 0.2, g_strand["q1"], g_strand["q2"], line_color="black")

    p.rect(["C Strand"], c_strand["lower"], whisker_width, c_strand["hbar_height"], line_color="grey")
    p.rect(["C Strand"], c_strand["upper"], whisker_width, c_strand["hbar_height"], line_color="grey")
    p.segment(["C Strand"], c_strand["upper"], ["C Strand"], c_strand["q3"], line_color="grey")
    p.segment(["C Strand"], c_strand["lower"], ["C Strand"], c_strand["q1"], line_color="grey")
    p.vbar(["C Strand"], 0.2, c_strand["q2"], c_strand["q3"], line_color="black")
    p.vbar(["C Strand"], 0.2, c_strand["q1"], c_strand["q2"], line_color="black")

    p.xaxis.major_label_orientation = "vertical"
    p.yaxis.axis_label = 'Values'


    if plt_title is not None:
        p.title.text = plt_title
    if x_title is not None:
        p.xaxis.axis_label = x_title
    if y_title is not None:
        p.yaxis.axis_label = y_title

    return plt




def telo_length_hist_by_strand(data_list=None, *, labels = None, x=None, y=None, hue=None, weights=None,
    stat='count', bins='auto', binwidth=None, binrange=None,
    discrete=None, cumulative=False, common_bins=True,
    common_norm=True, multiple='layer', element='bars',
    fill=True, shrink=1, kde=False, kde_kws=None,
    line_kws=None, thresh=0, pthresh=None, pmax=None,
    cbar=False, cbar_ax=None, cbar_kws=None, palette=None,
    hue_order=None, hue_norm=None, color=None, log_scale=None,
        legend=True, ax=None, plt_title=None, x_title=None, y_title=None, **quad_kwargs):
    """Plot univariate or multivariate histograms."""
    plt = BokehPlot()

    plot = figure()

    estimate_kws = dict(
        stat=stat,
        bins=bins,
        binwidth=binwidth,
        binrange=binrange,
        discrete=discrete,
        cumulative=cumulative,
    )

    #data = pd.DataFrame(data)

    estimator = Histogram(**estimate_kws)

        # multivariate data
    opacity = 0.5
    if palette is None:
        palette = util.choose_palette()

    # this just looks over values if data is 1D
    # for var, color in zip(data, cycle(palette)):
    for col, color, labels in zip(data_list, cycle(palette), labels):
        quad_kwargs = {}
        
        quad_kwargs["legend_label"] = labels
        variable_data = col.dropna()
        heights, edges = estimator(variable_data, weights=weights)
        new_heights = [height/sum(heights) * 100 for height in heights]
        plt._fig.quad(
            top=new_heights, bottom=0, left=edges[:-1], right=edges[1:],
            fill_color=color, fill_alpha=opacity, line_color=color, **quad_kwargs
        )
    
    #plt._fig.vspan(x=[np.mean(data)], line_width=[1], color="red")
    plt._fig.y_range.start = 0
    hover = plt._fig.select(dict(type=HoverTool))
    hover.tooltips = [(stat.capitalize(), "@top")]

    if plt_title is not None:
        plt._fig.title.text = plt_title
    if x_title is not None:
        plt._fig.xaxis.axis_label = x_title
    if y_title is not None:
        plt._fig.yaxis.axis_label = y_title

    return plt


def repeat_freq_histogram(data=None, *, x=None, y=None, z=None, hue=None, weights=None,
    stat='count', bins='auto', binwidth=None, binrange=None,
    discrete=None, cumulative=False, common_bins=True,
    common_norm=True, multiple='layer', element='bars',
    fill=True, shrink=1, kde=False, kde_kws=None,
    line_kws=None, thresh=0, pthresh=None, pmax=None,
    cbar=False, cbar_ax=None, cbar_kws=None, palette=None,
    hue_order=None, hue_norm=None, color=None, log_scale=None,
        legend=True, ax=None, plt_title=None, x_title=None, y_title=None, **quad_kwargs):
    """Plot univariate or multivariate histograms."""
    plt = BokehPlot()

    plot = figure()
    
    # print(plot)
    # print(type(plot))
    # glyph = VSpan(x=np.mean(data))
    # plot.add_glyph(glyph)

    # print(plt._fig)
    # print(type(plt._fig))

    estimate_kws = dict(
        stat=stat,
        bins=bins,
        binwidth=binwidth,
        binrange=binrange,
        discrete=discrete,
        cumulative=cumulative,
    )
    
    data = pd.DataFrame(data)

    estimator = Histogram(**estimate_kws)

    if data.ndim > 1 and data.shape[1] > 1:
        # multivariate data
        opacity = 0.5
        if palette is None:
            palette = util.choose_palette()
    else:
        opacity = 1.0
        if color is None:
            palette = util.choose_palette()
        else:
            palette = [color]
    if hue:
        data = data.pivot(columns=hue, values=data.columns[0])
    # this just looks over values if data is 1D
    # for var, color in zip(data, cycle(palette)):

    for col, color in zip(data.columns, cycle(palette)):
        # print(col)
        # print(color)
        quad_kwargs = {}
        if len(data.columns) > 1:
            quad_kwargs["legend_label"] = col
        variable_data = data[col].dropna()
        heights, edges = estimator(variable_data, weights=weights)

        plt._fig.quad(
            top=heights, bottom=0, left=edges[:-1], right=edges[1:],
            fill_color=color, fill_alpha=opacity, line_color=color, **quad_kwargs
        )
    
    #plt._fig.vspan(x=[np.mean(data)], line_width=[1], color="red")
    plt._fig.y_range.start = 0
    hover = plt._fig.select(dict(type=HoverTool))
    hover.tooltips = [(stat.capitalize(), "@top")]

    if plt_title is not None:
        plt._fig.title.text = plt_title
    if x_title is not None:
        plt._fig.xaxis.axis_label = x_title
    if y_title is not None:
        plt._fig.yaxis.axis_label = y_title

    return plt


def colored_telo_length_barplot(data=None, *, x=None, y=None, hue=None, order=None, hue_order=None,
    estimator='mean', errorbar=('ci', 95), n_boot=1000, units=None, seed=None,
    orient=None, color=None, palette=None, saturation=1.0, width=0.8,
    errcolor='.26', errwidth=None, capsize=None, dodge=True, ci='deprecated',
    ax=None, nested_x=False,  mutant=False, plt_title=None, y_title=None, x_title=None, x_rotation=None,
    repeat="GGTTAG", **kwargs,
):
    if not mutant:
        read_ids = data["read_id"]
        fill = ["Non-{} or One Nucl. Substitutions".format(repeat), "One Nucl. Substitition of {}".format(repeat), repeat]
        colors=["#D81B60", "#1E88E5","#FFC107"]
        source_data = {"read_ids" : read_ids,
                    "Non-{} or One Nucl. Substitutions".format(repeat): (100 - (data["one_nucl_variant_composition"] + data["wt_composition"]))/100 *data["vrr_telo_length"],
                    "One Nucl. Substitition of {}".format(repeat) : data["one_nucl_variant_composition"]/100*data["vrr_telo_length"],
                    repeat : data["wt_composition"]/100*data["vrr_telo_length"]}
        
        plt = BokehPlot(y_range=read_ids)
        p = plt._fig
        p.hbar_stack(fill, y="read_ids", height=0.9, color=colors, source=source_data, legend_label=fill)

    else:
        read_ids = data["read_id"]
     
        fill = ["Non-{} or One Nucl. Substitutions".format(repeat), "One Nucl. Substitition of {}".format(repeat), mutant, repeat]

        colors=["#D81B60", "#1E88E5","#FFC107", "#004D40"]
     
        source_data = {"read_ids" : read_ids,
                    "Non-{} or One Nucl. Substitutions".format(repeat): (100 - (data["one_nucl_variant_composition"] + data["wt_composition"]))/100 *data["vrr_telo_length"],
                    "One Nucl. Substitition of {}".format(repeat) : data["one_nucl_variant_composition"]/100*data["vrr_telo_length"],
                    repeat : data["wt_composition"]/100*data["vrr_telo_length"],
                    mutant : data["mutant_composition"]/100*data["vrr_telo_length"]}
        plt = BokehPlot(y_range=read_ids)
        p = plt._fig
        p.hbar_stack(fill, y="read_ids", height=0.9, color=colors, source=source_data, legend_label=fill)

    #p.yaxis.major_grid_line_color = None
    #p.yaxis.minor_grid_line_color = None
    p.ygrid.grid_line_color = None
    #p.yaxis.major_label_overrides([i for i in range(0, len(read_ids))])
    p.yaxis.visible = False
    p.legend.location = "bottom_right"
    if plt_title is not None:
        p.title.text = plt_title
    if x_title is not None:
        p.xaxis.axis_label = x_title
    if y_title is not None:
        p.yaxis.axis_label = y_title
    if x_rotation is not None:
        p.xaxis.major_label_orientation = x_rotation

    return plt