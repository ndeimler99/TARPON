#!/usr/bin/env python3
from collections import Counter

from bokeh.models import ColumnDataSource, FactorRange, Whisker
import bokeh.transform as bkt
import numpy as np
import pandas as pd
from seaborn.categorical import _BarPlotter

from ezcharts.plots import BokehPlot


def seqkit_stats_boxplot_length(
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
    iqr = df.Q3 - df.Q1

    # whiskers
    if isinstance(whis, float):
        df["upper"] = df.Q3 + whis * iqr
        df["lower"] = df.Q1 - whis * iqr
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
        base=x, upper="upper", lower="min_len", source=source,
        line_width=linewidth, line_color=linecolor)

    whisker.upper_head.size = whisker.lower_head.size = 20

    y_max = df["upper"].max()

    y_range = (0, 1.1*y_max)
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

    return plt