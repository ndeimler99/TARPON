#!/usr/bin/env python3

import sys
import os
import pandas as pd
#import plotly.express as px
#import plotly.graph_objects as go
import math
#from jinja2 import Template
import argparse
import json

#from bokeh.models import HoverTool
from dominate import tags as html_tags
from dominate.tags import em, p
import ezcharts as ezc
from ezcharts.components.ezchart import EZChart
from ezcharts.components.fastcat import SeqSummary
from ezcharts.components.reports.labs import LabsReport
from ezcharts.layout.snippets import DataTable
from ezcharts.layout.snippets import Grid
from ezcharts.layout.snippets import Tabs
import pandas as pd
from ezcharts.plots import BokehPlot
from ezcharts.components.reports.comp import ComponentReport
from bokeh.models import HoverTool
from report_utils import *


THEME = 'epi2melabs'

def main(args):

    params_file = open(args.params)
    params = params_file.read()
    params = json.loads(params)

    manifest_file = open(args.manifest)
    manifest = manifest_file.read()
    manifest = json.loads(manifest)

    versions = pd.read_table(args.versions, sep=",", header=None)
    versions.set_axis(["Software", "Version"], axis=1)

    sample_dict = {}
    for i in args.sample_stats_retained:
        sample_dict[i.split(".")[0]] = {}
        sample_dict[i.split(".")[0]]["retained_reads"] = i

    for i in args.sample_stats_filtered:
        sample_dict[i.split(".")[0]]["filtered_reads"] = i

    # versions_tables = go.Figure(data=[go.Table(
    #     header=dict(values=list(versions.columns),
    #             fill_color='paleturquoise',
    #             align='center'),
    #     cells=dict(values=versions.transpose().values.tolist(),
    #            fill_color='lavender',
    #            align='left'))
    # ])

    report = LabsReport(
        f"{args.workflow_name} Report for Run: {params['run_name']}", args.workflow_name,
        args.params, args.versions, manifest["version"])
    
    with report.add_section("Sequencing Stats", "Sequencing Stats"):
        p("""Statistics for the entire flow cell of non-demultiplexed data""")
        tabs = Tabs()
        with tabs.add_dropdown_menu("Retained Read Statistics", change_header=False):
            df = pd.read_table(args.run_stats_retained)
            df["file"] = df["file"].apply(lambda x: x.split("/")[-1].split(".fastq")[0])
            df = df.drop(columns=["N50.1", "format", "sum_len", "sum_gap"])
            with tabs.add_dropdown_tab('Stats Table'):
                DataTable.from_pandas(df, use_index=False)
            with tabs.add_dropdown_tab('Read Count'):
                plt = ezc.barplot(data=df, x="file", y="num_seqs")
                plt._fig.title.text = "Number of Reads Retained at Each Step"
                plt._fig.xaxis.axis_label = "Pipeline Step"
                plt._fig.yaxis.axis_label = "Number of Retained Reads"
                plt._fig.xaxis.major_label_orientation = 45
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Number of Reads Retained", "@top")]
                print(plt)
                EZChart(plt,THEME)

            with tabs.add_dropdown_tab('Read Length'):
                plt = seqkit_stats_boxplot_length(df, x="file")
                plt._fig.yaxis.axis_label = "Read Length (BP)"
                plt._fig.xaxis.axis_label = "Pipeline Step"
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Length", "@avg_len"),
                                  ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                   ("Min Length", "@min_len"), ("Max Length", "@max_len")]
                EZChart(plt, THEME)

            with tabs.add_dropdown_tab('Read Quality'):
                plt=ezc.barplot(data=df, x="file", y="AvgQual")
                plt._fig.title.text = "Average Quality of Reads Retained at Each Step"
                plt._fig.xaxis.axis_label = "Pipeline Step"
                plt._fig.yaxis.axis_label = "Average Quality"
                plt._fig.xaxis.major_label_orientation = 45
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Average Quality", "@top"), ("Read Count", "@num_seqs")]
                EZChart(plt,THEME)

        with tabs.add_dropdown_menu("Filtered Read Statistics", change_header=False):
            df = pd.read_table(args.run_stats_filtered)
            df["file"] = df["file"].apply(lambda x: x.split("/")[-1].split(".fastq")[0])
            df = df.drop(columns=["N50.1", "format", "sum_len", "sum_gap"])

            with tabs.add_dropdown_tab('Stats Table'):
                DataTable.from_pandas(df, use_index=False)

            with tabs.add_dropdown_tab('Read Count'):
                plt=ezc.barplot(data=df, x="file", y="num_seqs")
                plt._fig.title.text = "Number of Reads Retained at Each Step"
                plt._fig.xaxis.axis_label = "Pipeline Step"
                plt._fig.yaxis.axis_label = "Number of Filtered Reads"
                plt._fig.xaxis.major_label_orientation = 45
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Number of Reads Filtered", "@top")]
                EZChart(plt,THEME)
            with tabs.add_dropdown_tab('Read Length'):
                plt = seqkit_stats_boxplot_length(df, x="file")
                plt._fig.yaxis.axis_label = "Read Length (BP)"
                plt._fig.xaxis.axis_label = "Pipeline Step"
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Length", "@avg_len"),
                                  ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                   ("Min Length", "@min_len"), ("Max Length", "@max_len")]
                EZChart(plt, THEME)

            with tabs.add_dropdown_tab('Read Quality'):
                plt=ezc.barplot(data=df, x="file", y="AvgQual")
                plt._fig.title.text = "Average Quality of Reads Filtered at Each Step"
                plt._fig.xaxis.axis_label = "Pipeline Step"
                plt._fig.yaxis.axis_label = "Average Quality"
                plt._fig.xaxis.major_label_orientation = 45
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Average Quality", "@top"), ("Read Count", "@num_seqs")]
                EZChart(plt,THEME)

    with report.add_section("Sample Comparison", "Sample Comparison"):
        p("""Between Sample Comparisons of Telomere Length""")
        tabs = Tabs()
        with tabs.add_tab("Number of Telomeric Reads"):
            pass
        with tabs.add_tab("Telo Stats"):
            pass
        with tabs.add_tab("Telomere Length Barchart"):
            pass
        with tabs.add_tab("Telomere Length Boxplot"):
            pass

    with report.add_section("Individual Sample Analysis", "Individual Sample Analysis"):
        p("Individual Sample Statistics and Plots")
        tabs = Tabs()
        for sample in sample_dict:
            with tabs.add_dropdown_menu(sample, change_header=False):
                with tabs.add_dropdown_tab("Retained Reads"):
                    pass
                with tabs.add_dropdown_tab("Filtered Reads"):
                    pass
                with tabs.add_dropdown_tab("Telomere Length"):
                    pass
                #if restriction digest
                #if vrr
                #if detailed
                


    # this works for optional/additional analysis
    # if args.restriction_digest_analysis is not None:
    #     with report.add_section("Restriction Digest Analysis", "Restriction Digest Analysis"):
    #         pass

    report.write(args.report)

def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", required=True)
    parser.add_argument("--template_file", required=True)
    parser.add_argument("--workflow_name", required=True, help="The name of the workflow.") # works
    parser.add_argument("--params", required=True) #works
    parser.add_argument("--versions", required=True) #works  
    parser.add_argument("--manifest", required=True) #works
    parser.add_argument("--commandLine", required=True)
    parser.add_argument("--run_stats_retained", required=True) #works but need to figure out how to get combined sample stats here
    parser.add_argument("--run_stats_filtered", required=True) #works but need to figure out how to get combined sample stats here
    parser.add_argument("--sample_stats_retained", nargs='+', required=True) #works but only for simplex, not multiplex tested yet
    parser.add_argument("--sample_stats_filtered", nargs='+', required=True) #works but only for simplex, not multiplex tested yet
    parser.add_argument("--restriction_digest_analysis", required=False, default="test")
    
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)

