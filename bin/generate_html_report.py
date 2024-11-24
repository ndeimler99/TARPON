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
from report_utils import BokehPlot
from ezcharts.components.reports.comp import ComponentReport
from bokeh.models import HoverTool, ColumnDataSource
import report_utils
import numpy as np
from bokeh.plotting import figure
import bokeh.palettes


# reorder statistics so input is first...
#  Individual Sample histogram and boxplot add red line denoting mean
#  Individual sample boxplot and barchart remove x axis tick labels
#  Add Manifest tablel similar to params table
#  Report page – remove provided by Oxford Nanopore Technologies
#  Report page – remove About footer – or at least modify text
#  Report page – remove Epi2Me logo and replace with TArPON cartoon
#  change read quality to boxplot



THEME = 'epi2melabs'

def get_nextflow_attributes(attribute_file):

    attribute_file = open(attribute_file)
    attribute = attribute_file.read()
    attribute = json.loads(attribute)
    return attribute

def main(args):

    args.plot_telo_length = args.plot_telo_length == "true"
    args.plot_vrr_length = args.plot_vrr_length == "true"
    args.strand_comparison = args.strand_comparison == "true"
    args.detailed_stats = args.detailed_stats == "true"

    params = get_nextflow_attributes(args.params)
    manifest = get_nextflow_attributes(args.manifest)

    versions = pd.read_table(args.versions, sep=",", header=None)
    versions.set_axis(["Software", "Version"], axis=1)

    sample_dict = {}
    for i in args.sample_stats_retained:
        sample_dict[i.split(".")[0]] = {}
        sample_dict[i.split(".")[0]]["retained_reads"] = i

    for i in args.sample_stats_filtered:
        sample_dict[i.split(".")[0]]["filtered_reads"] = i

    for i in args.sample_telo_stats:
        sample_dict[i.split(".")[0]]["telo_stats"] = i

    if args.restriction_digest[0] != "false":
        if i.endswith(".txt"):
            for i in args.restriction_digest:
                sample_dict[i.split(".")[0]]["digest"] = i
    

    if args.plot_telo_length:
        telo_summary_stats = pd.read_table(args.run_telo_stats, sep="\t")

    if args.plot_vrr_length:
        vrr_summary_stats = pd.read_table(args.run_vrr_stats, sep="\t")

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
            with tabs.add_dropdown_tab('Read Count'):
                plt = report_utils.barplot(data=df, x="file", y="num_seqs", 
                                           x_title="Pipeline Step", x_rotation=45,
                                           y_title="Number of Retained_Reads",
                                           plt_title="Number of Reads Retained at Each Step")
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Number of Reads Retained", "@num_seqs")]
                # what is this function doing?
                EZChart(plt,THEME)
            with tabs.add_dropdown_tab('Stats Table'):
                DataTable.from_pandas(df, use_index=False)
            with tabs.add_dropdown_tab('Read Length'):
                plt = report_utils.seqkit_stats_boxplot_length(df, x="file", 
                                                               y_title="Read Length (BP)", x_title="Pipeline Step", 
                                                               x_rotation=45,
                                                               plt_title="Read Length of Reads Retained at Each Step")
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Length", "@avg_len"),
                                  ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                   ("Min Length", "@min_len"), ("Max Length", "@max_len")]
                # what is this function doing?
                EZChart(plt, THEME)

            with tabs.add_dropdown_tab('Read Quality'):
                plt = report_utils.barplot(data=df, x="file", y="AvgQual",
                                           x_title="Pipeline Step", y_title="Average Quality",
                                           plt_title="Average Quality of Reads Retained at Each Step",
                                           x_rotation=45)
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("File", "@file"), ("Average Quality", "@AvgQual"), ("Read Count", "@num_seqs")]
                EZChart(plt,THEME)

        with tabs.add_dropdown_menu("Filtered Read Statistics", change_header=False):
            df = pd.read_table(args.run_stats_filtered)
            df["file"] = df["file"].apply(lambda x: x.split("/")[-1].split(".fastq")[0])
            df = df.drop(columns=["N50.1", "format", "sum_len", "sum_gap"])

            with tabs.add_dropdown_tab('Read Count'):
                plt=report_utils.barplot(data=df, x="file", y="num_seqs",
                                plt_title="Number of Reads Filtered at Each Step",
                                x_title="Pipeline Step", x_rotation=45,
                                y_title="Number of Filtered Reads")
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Number of Reads Filtered", "@num_seqs")]
                # what is this function doing
                EZChart(plt,THEME)
            with tabs.add_dropdown_tab('Stats Table'):
                DataTable.from_pandas(df, use_index=False)
            with tabs.add_dropdown_tab('Read Length'):
                plt = report_utils.seqkit_stats_boxplot_length(df, x="file",
                                                               y_title="Read Length (BP)", 
                                                               x_title="Pipeline Step", x_rotation=45,
                                                               plt_title="Read Length of Reads Filtered at Each Step")
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Length", "@avg_len"),
                                  ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                   ("Min Length", "@min_len"), ("Max Length", "@max_len")]
                # what is this function doing?
                EZChart(plt, THEME)

            with tabs.add_dropdown_tab('Read Quality'):
                plt=report_utils.barplot(data=df, x="file", y="AvgQual",
                                         x_title="Pipeline Step", y_title="Average Quality",
                                         plt_title="Average Quality of Reads Filtered at Each Step",
                                         x_rotation=45)
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("File", "@file"), ("Average Quality", "@AvgQual"), ("Read Count", "@num_seqs")]
                # what is this function doing
                EZChart(plt,THEME)

    with report.add_section("Sample Comparison", "Sample Comparison"):
        p("""Between Sample Comparisons of Telomere Length""")
        tabs = Tabs()
        with tabs.add_tab("Number of Telomeric Reads"):
            #barplot using numseqs
            if args.plot_telo_length:
                df = telo_summary_stats
            else:
                df = vrr_summary_stats

            df = df.set_index("Sample_ID").loc[sorted(sample_dict.keys())].reset_index()
            plt=report_utils.barplot(data=df, x="Sample_ID", y="Number_of_Reads",
                                     plt_title="Number of Telomeric Reads per Sample",
                                     x_title="Sample", x_rotation=45,
                                     y_title="Number of Telomeres")
            hover = plt._fig.select(dict(type=HoverTool))
            if args.plot_telo_length:
                hover.tooltips = [("Sample", "@Sample_ID"), ("Number of Reads", "@Number_of_Reads"),("Mean_Telomere_Length", "@Mean_Telomere_Length")]
            else:
                hover.tooltips = [("Sample", "@Sample_ID"), ("Number of Reads", "@Number_of_Reads"),("Mean_Telomere_Length", "@Mean_VRR_Telomere_Length")]

            EZChart(plt, THEME)
            ####
        if args.plot_telo_length:
            with tabs.add_tab("Telo Stats"):
                df = telo_summary_stats
                df = df.set_index("Sample_ID").loc[sorted(sample_dict.keys())].reset_index()
                DataTable.from_pandas(df, use_index=False)
            with tabs.add_tab("Telomere Length Barchart"):
                master_df = pd.DataFrame()
                for sample in sample_dict.keys():
                    df = pd.read_table(sample_dict[sample]["telo_stats"], sep="\t")
                    bins = [i*1000 for i in range(0,11)]
                    bins.append(100000)
                    telo_bar_df = np.histogram(df["telo_length"], bins=bins)
                    telo_bar_df = pd.DataFrame(list(zip(telo_bar_df[1], telo_bar_df[0])), columns=["bin_start", "bin_size"])
                    telo_bar_df["bin_start"] = telo_bar_df["bin_start"].astype("string")
                    telo_bar_df["sample"] = sample
                    telo_bar_df["bin_size"] = telo_bar_df["bin_size"] / sum(telo_bar_df["bin_size"]) * 100
                    master_df = pd.concat([master_df, telo_bar_df], ignore_index=True, sort=False)
                

                telo_bar_plot = report_utils.telo_barplot(data=master_df, 
                                                    x="sample", x_rotation=45, x_title="Sample",
                                                    y="bin_size", y_title="Percentage of Reads",
                                                    hue="bin_start", dodge=False,
                                                    order=sorted(list(sample_dict.keys())), 
                                                    palette=bokeh.palettes.Category20[11],
                                                    plt_title="Telomere Length by Sample Binned Bar Plot")
                EZChart(telo_bar_plot, THEME)
            with tabs.add_tab("Telomere Length Boxplot"):
                df = telo_summary_stats
                plt = report_utils.seqkit_stats_boxplot_length(df, x="Sample_ID",
                                                            plt_title="Telomere Length Boxplots",
                                                            x_title="Sample", x_rotation=45,
                                                            y_title="Telomere Length (BP)",
                                                            order=sorted(list(sample_dict.keys())))
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Sample", "@Sample_ID"),("Number of Reads", "@Number_of_Reads"),
                                ("Avg Length", "@Mean_Telomere_Length"),
                                    ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                    ("Min Length", "@Min_Telo_Length"), ("Max Length", "@Max_Telo_Length")]
                EZChart(plt,THEME)
        
        if args.plot_vrr_length:
            with tabs.add_tab("VRR Stats"):
                df = vrr_summary_stats
                df = df.set_index("Sample_ID").loc[sorted(sample_dict.keys())].reset_index()
                DataTable.from_pandas(df, use_index=False)
            with tabs.add_tab("VRR Length Barchart"):
                master_df = pd.DataFrame()
                for sample in sample_dict.keys():
                    df = pd.read_table(sample_dict[sample]["telo_stats"], sep="\t")
                    bins = [i*1000 for i in range(0,11)]
                    bins.append(100000)
                    telo_bar_df = np.histogram(df["vrr_telo_length"], bins=bins)
                    telo_bar_df = pd.DataFrame(list(zip(telo_bar_df[1], telo_bar_df[0])), columns=["bin_start", "bin_size"])
                    telo_bar_df["bin_start"] = telo_bar_df["bin_start"].astype("string")
                    telo_bar_df["sample"] = sample
                    telo_bar_df["bin_size"] = telo_bar_df["bin_size"] / sum(telo_bar_df["bin_size"]) * 100
                    master_df = pd.concat([master_df, telo_bar_df], ignore_index=True, sort=False)
                

                telo_bar_plot = report_utils.telo_barplot(data=master_df, 
                                                    x="sample", x_rotation=45, x_title="Sample",
                                                    y="bin_size", y_title="Percentage of Reads",
                                                    hue="bin_start", dodge=False,
                                                    order=sorted(list(sample_dict.keys())), 
                                                    palette=bokeh.palettes.Category20[11],
                                                    plt_title="VRR Telomere Length by Sample Binned Bar Plot")
                EZChart(telo_bar_plot, THEME)
            with tabs.add_tab("VRR Length Boxplot"):
                df = vrr_summary_stats
                plt = report_utils.seqkit_stats_boxplot_length(df, x="Sample_ID",
                                                            plt_title="VRR Telomere Length Boxplots",
                                                            x_title="Sample", x_rotation=45,
                                                            y_title="VRR Telomere Length (BP)",
                                                            order=sorted(list(sample_dict.keys())))
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Sample", "@Sample_ID"),("Number of Reads", "@Number_of_Reads"),
                                ("Avg Length", "@Mean_VRR_Telomere_Length"),
                                    ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                    ("Min Length", "@Min_VRR_Telo_Length"), ("Max Length", "@Max_VRR_Telo_Length")]
                EZChart(plt,THEME)

    with report.add_section("Individual Sample Analysis", "Individual Sample Analysis"):
        # p("Individual Sample Statistics and Plots")
        tabs = Tabs()
        for sample in sorted(list(sample_dict.keys())):
            with tabs.add_dropdown_menu(sample, change_header=False):
                with tabs.add_dropdown_tab("{} Telomere Length".format(sample)):
                    df = pd.read_table(sample_dict[sample]["telo_stats"], sep="\t")
                    df["telo_length"] = df["telo_length"].astype("float")
                    #telo length histogram next to barplot
                    new_tabs = Tabs()
                    if args.plot_telo_length:
                        with new_tabs.add_tab("Telo Length Analysis (n={})".format(len(df['telo_length']))):
                            # Telomere Length Histogram
                            telo_length_hist = report_utils.telo_length_hist(df["telo_length"], binwidth=200, binrange=[0,max(df["telo_length"])+200],
                                                            plt_title="Telomere Length Histogram",
                                                            x_title = "Telomere Length (BP)", y_title="Read Count")
                        
                            # Telomere Length Boxplot
                            telo_length_boxplot = report_utils.create_boxplot(df=df,column_name="telo_length", sample=sample,
                                                                            x_title = sample, y_title="Telomere Length",
                                                                            plt_title="Telomere Length Boxplot")
                
                            # Telomere Length Barplot
                            bins = [i*1000 for i in range(0,11)]
                            bins.append(100000)
                            telo_bar_df = np.histogram(df["telo_length"], bins=bins)
                            telo_bar_df = pd.DataFrame(list(zip(telo_bar_df[1], telo_bar_df[0])), columns=["bin_start", "bin_size"])
                            telo_bar_df["bin_start"] = telo_bar_df["bin_start"].astype("string")
                            telo_bar_df["sample"] = sample
                            telo_bar_df["bin_size"] = telo_bar_df["bin_size"] / sum(telo_bar_df["bin_size"]) * 100
                            #telo_bar_df["bin_label"] = "{} - {} bp".format(telo_bar_df["bin_start"] + 1, telo_bar_df["bin_start"]+1000)
                            telo_bar_plot = report_utils.telo_barplot(data=telo_bar_df, x="sample", 
                                                                    y="bin_size", hue="bin_start", dodge=False, 
                                                                    palette=bokeh.palettes.Category20[11],
                                                                    plt_title="Telomere Length Binned",
                                                                    x_title=sample, y_title="Percentage of Telomeres", 
                                                                    legend_loc="right", legend_orientation="vertical",
                                                                    hide_x_tick_labels=True)
                            with Grid(columns=3):
                                EZChart(telo_length_hist, THEME)
                                EZChart(telo_length_boxplot, THEME)
                                EZChart(telo_bar_plot, THEME)

                        with new_tabs.add_tab("Telo Length vs Read Length"):
                            telo_length_hist = report_utils.telo_length_hist(df["telo_length"], binwidth=200, binrange=[0,max(df["telo_length"])+200],
                                                                            plt_title="Telomere Length Histogram",
                                                                            x_title="Telomere Length (BP)", y_title="Read Count")
        
                            read_length_hist = report_utils.telo_length_hist(df["read_len"], binwidth=500, binrange=[0, max(df["read_len"])+500],
                                                                            plt_title="Read Length Histogram", x_title="Read Length (BP)",
                                                                            y_title="Read Count")
    
                            byscatter = report_utils.scatterplot(data=df, x="read_len", y="telo_length", 
                                                                hover_tooltips=[("Telomere Length", "@y"), ("Read Length", "@x")],
                                                                plt_title="Telomere Length by Read Length",
                                                                x_title="Read Length (BP)", y_title="Telomere Length (BP)")
                            # byscatter.xAxis.name = "Read Length"
                            # byscatter.yAxis.name = "Telomere Length"
                            with Grid(columns=3):
                                #telo length hist, read length_hist, scatter plot
                                EZChart(telo_length_hist, THEME)
                                EZChart(read_length_hist, THEME)
                                EZChart(byscatter, THEME)
                    if args.plot_vrr_length:
                        with new_tabs.add_tab("VRR Length Analysis (n={})".format(len(df['telo_length']))):
                            # Telomere Length Histogram
                            telo_length_hist = report_utils.telo_length_hist(df["vrr_telo_length"], binwidth=200, binrange=[0,max(df["vrr_telo_length"])+200],
                                                            plt_title="VRR Telomere Length Histogram",
                                                            x_title = "VRR Telomere Length (BP)", y_title="Read Count")
                        
                            # Telomere Length Boxplot
                            telo_length_boxplot = report_utils.create_boxplot(df=df,column_name="vrr_telo_length", sample=sample,
                                                                            x_title = sample, y_title="VRR Telomere Length",
                                                                            plt_title="VRR Telomere Length Boxplot")
                
                            # Telomere Length Barplot
                            bins = [i*1000 for i in range(0,11)]
                            bins.append(100000)
                            telo_bar_df = np.histogram(df["vrr_telo_length"], bins=bins)
                            telo_bar_df = pd.DataFrame(list(zip(telo_bar_df[1], telo_bar_df[0])), columns=["bin_start", "bin_size"])
                            telo_bar_df["bin_start"] = telo_bar_df["bin_start"].astype("string")
                            telo_bar_df["sample"] = sample
                            telo_bar_df["bin_size"] = telo_bar_df["bin_size"] / sum(telo_bar_df["bin_size"]) * 100
                            #telo_bar_df["bin_label"] = "{} - {} bp".format(telo_bar_df["bin_start"] + 1, telo_bar_df["bin_start"]+1000)
                            telo_bar_plot = report_utils.telo_barplot(data=telo_bar_df, x="sample", 
                                                                    y="bin_size", hue="bin_start", dodge=False, 
                                                                    palette=bokeh.palettes.Category20[11],
                                                                    plt_title="VRR Telomere Length Binned",
                                                                    x_title=sample, y_title="Percentage of Telomeres", 
                                                                    legend_loc="right", legend_orientation="vertical",
                                                                    hide_x_tick_labels=True)
                            with Grid(columns=3):
                                EZChart(telo_length_hist, THEME)
                                EZChart(telo_length_boxplot, THEME)
                                EZChart(telo_bar_plot, THEME)
                        with new_tabs.add_tab("VRR Length vs Read Length"):
                            telo_length_hist = report_utils.telo_length_hist(df["vrr_telo_length"], binwidth=200, binrange=[0,max(df["vrr_telo_length"])+200],
                                                                            plt_title="VRR Telomere Length Histogram",
                                                                            x_title="VRR Telomere Length (BP)", y_title="Read Count")
        
                            read_length_hist = report_utils.telo_length_hist(df["read_len"], binwidth=500, binrange=[0, max(df["read_len"])+500],
                                                                            plt_title="Read Length Histogram", x_title="Read Length (BP)",
                                                                            y_title="Read Count")
    
                            byscatter = report_utils.scatterplot(data=df, x="read_len", y="vrr_telo_length", 
                                                                hover_tooltips=[("VRR Telomere Length", "@y"), ("Read Length", "@x")],
                                                                plt_title="VRR Telomere Length by Read Length",
                                                                x_title="Read Length (BP)", y_title="VRR Telomere Length (BP)")
                            # byscatter.xAxis.name = "Read Length"
                            # byscatter.yAxis.name = "Telomere Length"
                            with Grid(columns=3):
                                #telo length hist, read length_hist, scatter plot
                                EZChart(telo_length_hist, THEME)
                                EZChart(read_length_hist, THEME)
                                EZChart(byscatter, THEME)
                    if args.plot_vrr_length and args.plot_telo_length:
                        with new_tabs.add_tab("VRR Length vs Telo Length"):
                            vrr_hist = report_utils.telo_length_hist(df["vrr_telo_length"], binwidth=200, binrange=[0, max(df["vrr_telo_length"])+200],
                                                                    plt_title="VRR Telomere Length",
                                                                    x_title="VRR Telomere Length (BP)", y_title="Read Count")
                            telo_hist = report_utils.telo_length_hist(df["telo_length"], binwidth=200, binrange=[0, max(df["telo_length"])+200],
                                                                    plt_title="Telomere LEngth", x_title = "Telomere Length (BP)",
                                                                    y_title="Read Count")
                            telo_scatter = report_utils.scatterplot(data=df, x="vrr_telo_length", y="telo_length",
                                                                    plt_title = "VRR Telomere Length By Telomere Length",
                                                                    x_title="VRR Telomere Length (BP)", y_title="Telomere Length (BP)",
                                                                    hover_tooltips=[("VRR Telomere Length", "@x"), ("Telomere Length", "@y")])
                            with Grid(columns=3):
                                EZChart(vrr_hist, THEME)
                                EZChart(telo_hist, THEME)
                                EZChart(telo_scatter, THEME)                                    
                with tabs.add_dropdown_tab("{} Retained Reads".format(sample)):
                    new_tabs = Tabs()
                    df = pd.read_table(sample_dict[sample]["retained_reads"], sep="\t")
                    df["file"] = df["file"].apply(lambda x: x.split("/")[-1].split(".fastq")[0])
                    df = df.drop(columns=["N50.1", "format", "sum_len", "sum_gap"])
                    with new_tabs.add_tab("Read Length"):
                        #bar plot
                        plt = report_utils.barplot(data=df, x="file", y="num_seqs",
                                                   plt_title="Number of Reads Retained at Each Step",
                                                   x_title="Pipeline Step", x_rotation=45, y_title="Number of Retained Reads")
                        hover = plt._fig.select(dict(type=HoverTool))
                        hover.tooltips = [("Number of Reads Retained", "@num_seqs")]
                        EZChart(plt,THEME)
                    #stats
                    with new_tabs.add_tab("Read Stats"):
                        DataTable.from_pandas(df, use_index=False)
                    #length boxplot
                    with new_tabs.add_tab("Read Length"):
                        plt = report_utils.seqkit_stats_boxplot_length(df, x="file",
                                                                       plt_title="Read Length of Retained Reads",
                                                                       x_title="Pipeline Step", x_rotation=45,
                                                                       y_title="Read Length (BP)")
                        hover = plt._fig.select(dict(type=HoverTool))
                        hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                        ("Avg Length", "@avg_len"),
                                        ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                        ("Min Length", "@min_len"), ("Max Length", "@max_len")]
                        EZChart(plt, THEME)
                    #quality plot
                    with new_tabs.add_tab("Read Quality"):
                        plt=report_utils.barplot(data=df, x="file", y="AvgQual",
                                                 plt_title="Average Quality of Reads Retained",
                                                 y_title="Average Quality",
                                                 x_title="Pipeline Step", x_rotation=45)
                        hover = plt._fig.select(dict(type=HoverTool))
                        hover.tooltips = [("File", "@file"), ("Average Quality", "@AvgQual"), ("Read Count", "@num_seqs")]
                        EZChart(plt,THEME)
                with tabs.add_dropdown_tab("{} Filtered Reads".format(sample)):
                    new_tabs = Tabs()
                    df = pd.read_table(sample_dict[sample]["filtered_reads"], sep="\t")
                    df["file"] = df["file"].apply(lambda x: x.split("/")[-1].split(".fastq")[0])
                    df = df.drop(columns=["N50.1", "format", "sum_len", "sum_gap"])
                    with new_tabs.add_tab("Read Length"):
                        #bar plot
                        plt = report_utils.barplot(data=df, x="file", y="num_seqs", 
                                                   plt_title="Number of Reads Filtered at Each Step",
                                                   x_title="Pipeline Step", x_rotation=45,
                                                   y_title="Number of Retained Reads")
                        hover = plt._fig.select(dict(type=HoverTool))
                        hover.tooltips = [("Number of Reads Filtered", "@num_seqs")]
                        EZChart(plt,THEME)
                    #stats
                    with new_tabs.add_tab("Read Stats"):
                        DataTable.from_pandas(df, use_index=False)
                    #length boxplot
                    with new_tabs.add_tab("Read Length"):
                        plt = report_utils.seqkit_stats_boxplot_length(df, x="file",
                                                                       y_title="Read Length (BP)",
                                                                       x_title="Pipeline Step",
                                                                       x_rotation=45)
                        hover = plt._fig.select(dict(type=HoverTool))
                        hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                        ("Avg Length", "@avg_len"),
                                        ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                        ("Min Length", "@min_len"), ("Max Length", "@max_len")]
                        EZChart(plt, THEME)
                    #quality plot
                    with new_tabs.add_tab("Read Quality"):
                        plt=report_utils.barplot(data=df, x="file", y="AvgQual", 
                                                 plt_title="Average Quality of Reads Filtered at Each Step",
                                                 y_title="Average Quality", x_title="Pipeline Step",
                                                 x_rotation=45)
                        hover = plt._fig.select(dict(type=HoverTool))
                        hover.tooltips = [("File", "@file"), ("Average Quality", "@AvgQual"), ("Read Count", "@num_seqs")]
                        EZChart(plt,THEME)
                if args.restriction_digest[0] != "false":
                        with tabs.add_dropdown_tab("{} Restriction Digest".format(sample)):
                            stats_df = pd.read_table(sample_dict[sample]["digest"], sep="\t")
                            stats_df = stats_df.drop(columns=["N50.1", "format", "sum_len", "sum_gap"])
                            DataTable.from_pandas(stats_df, use_index=False)

                if args.strand_comparison:
                
                # where to add some sort of descriptive stats? with actual numbers, not just figures
                    with tabs.add_dropdown_tab("Strand Comparison"):
                        new_tabs = Tabs()
                        with new_tabs.add_tab("Read Length"):
                            with Grid(columns=3):
                                # can I add an additional Grid to split one of the grids into two?
                #             # bar plot on number of reads
                #             # read length by strand
                #                 # G and C strand histograms one top of each other
                #             # boxplot of read length by strand
                                pass
                        if args.plot_telo_length:
                            with new_tabs.add_tab("Telo Length by Strand"):
                #                 #sample plots as read length but with telo length
                #                 # telo length histogram bar chart plot
                                pass
                        if args.plot_vrr_length:
                            with new_tabs.add_tab("VRR Telo Length by Strand"):
                #                 # see previous
                                pass
    
    if args.detailed_stats:
        with report.add_section("Detailed Analysis", "Detailed Analysis"):
            tabs = Tabs()
            for sample in sorted(list(sample_dict.keys())):
                print(sample_dict[sample]["telo_stats"])
                df = pd.read_table(sample_dict[sample]["telo_stats"], sep="\t")
                with tabs.add_dropdown_menu(sample, change_header=False):
                    with tabs.add_dropdown_tab("Sequencing Quality"):
                        new_tabs = Tabs()
                        with new_tabs.add_tab("Read vs. Telo Quality Comparison"):
                            # read quality hist
                            read_qual = report_utils.telo_length_hist(df["read_quality"], binwidth=1, 
                                                                    binrange=[0, max(df["read_quality"] + 1)],
                                                                    plt_title="Read Quality",
                                                                    x_title="Quality Score",
                                                                    y_title="Read Count")
            
                            # telo quality hist
                            telo_qual = report_utils.telo_length_hist(df["telo_quality"], binwidth=1, 
                                                                    binrange=[0, max(df["telo_quality"] + 1)],
                                                                    plt_title="Telomere Quality",
                                                                    x_title="Quality Score",
                                                                    y_title="Telomere Count")
                            # read quality vs telo scatterplot
                            read_telo_scatter = report_utils.scatterplot(data=df, x="telo_quality", y="read_quality",
                                                                        plt_title="Telomere vs Read Sequencing Quality",
                                                                        x_title= "Telomere Quality", y_title="Read Quality",
                                                                       hover_tooltips=[("Telomere Quality", "@x"), ("Read Quality", "@y")])
                            
                            
                            boxplot = report_utils.create_multisample_boxplot(df=df, column_names=["telo_quality", "read_quality"], 
                                                                                min_q = 15, max_q = 50, plt_title="Seq. Quality",
                                                                                y_title="Quality Score", x_rotation=45,
                                                                                x_labels={"telo_quality":"Telomere Quality", "read_quality":"Read Quality"})
                            
                            with Grid(columns=4):
                                EZChart(read_qual, THEME)
                                EZChart(telo_qual, THEME)            
                                EZChart(boxplot, THEME)
                                EZChart(read_telo_scatter, THEME)
                                
                        if args.plot_telo_length:
                            with new_tabs.add_tab("Quality vs Telomere Length"):
                                read_by_telo = report_utils.scatterplot(data=df, x="read_quality", y="telo_length",
                                                                        plt_title="Read Quality by Telo Length",
                                                                        x_title="Read Quality", y_title="Telomere Length (BP)",
                                                                        hover_tooltips=[("Read Quality", "@x"), ("Telomere Length", "@y")])
                                # read quality by telo length
                                # telo quality by telo length
                                telo_by_telo = report_utils.scatterplot(data=df, x="telo_quality", y="telo_length",
                                                                        plt_title="Telomere Quality by Telo Length",
                                                                        x_title="Telomere Quality", y_title="Telomere Length (BP)",
                                                                        hover_tooltips=[("Telomere Quality", "@x"), ("Telomere Length", "@y")])
                                
                                with Grid(columns=2):
                                    EZChart(read_by_telo, THEME)
                                    EZChart(telo_by_telo, THEME)
                        
                        if args.plot_vrr_length:
                            with new_tabs.add_tab("Quality vs VRR Length"):
                                # read quality by vrr length
                                read_by_telo = report_utils.scatterplot(data=df, x="read_quality", y="vrr_telo_length",
                                                                        plt_title="Read Quality by VRR Telo Length",
                                                                        x_title="Read Quality", y_title="VRR Telomere Length (BP)",
                                                                        hover_tooltips=[("Read Quality", "@x"), ("VRR Telomere Length", "@y")])
                                # telo quality by telo length
                                telo_by_telo = report_utils.scatterplot(data=df, x="telo_quality", y="vrr_telo_length",
                                                                        plt_title="Telomere Quality by VRR Telo Length",
                                                                        x_title="Telomere Quality", y_title="VRR Telomere Length (BP)",
                                                                        hover_tooltips=[("Telomere Quality", "@x"), ("VRR Telomere Length", "@y")])
                                
                                with Grid(columns=2):
                                    EZChart(read_by_telo, THEME)
                                    EZChart(telo_by_telo, THEME)
                        if args.strand_comparison:
                            with new_tabs.add_tab("Strand Comparison"):
                                # telomere perc boxplot comparison
                                # one sub perc boxplot comparison
                                # telo quality differences
                                # read quality boxplot differences
                                pass
                    with tabs.add_dropdown_tab("Telomeric Composition"):
                        new_tabs = Tabs()
                        with new_tabs.add_tab("Perfect Repeats Distribution"):
                                # telo % GGTTAG hist
                                perfect_hist = report_utils.telo_length_hist(df["perc_perfect"], binwidth=1, 
                                                                    binrange=[0, 100],
                                                                    plt_title="Telomeric Perfect Repeat Composition",
                                                                    x_title="Perentage of Perfect Repeats",
                                                                    y_title="Telomere Count")
                                # telo % GGTTAG boxplot
                                perfect_box = report_utils.create_boxplot(df=df,column_name="perc_perfect", sample=sample,
                                                                            x_title = sample, y_title="Percentage of Telomere",
                                                                            plt_title="Telomeric Perfect Repeat Composition")
                
                                # % GGTTAG by telomere quality
                                perc_by_quality = report_utils.scatterplot(data=df, x="telo_quality", y="perc_perfect",
                                                                            plt_title="Telomere Quality by Composition",
                                                                            x_title="Telomere Quality", 
                                                                            y_title="Percentage of Perfect Repeats",
                                                                            hover_tooltips=[("Telomere Quality", "@x"), ("Percentage Perfect Repeats", "@y")])
                                with Grid(columns=3):
                                    EZChart(perfect_hist, THEME)
                                    EZChart(perfect_box, THEME)
                                    EZChart(perc_by_quality, THEME)
                        with new_tabs.add_tab("Perfect Repeats by Telomere Length"):
                            # % GGTTAG by telomere length
                            length_by_perfect = report_utils.scatterplot(data=df, x="telo_length", y="perc_perfect",
                                                                            plt_title="Percentage of Perfect Repeats by Telomere Length",
                                                                            x_title="Telomere Length (BP)", y_title="Percentage of Perfect Repeats",
                                                                            hover_tooltips=[("Telomere Length", "@x"), ("Percentage Perfect Repeats", "@y")])
                            vrr_length_by_perfect = report_utils.scatterplot(data=df, x="vrr_telo_length", y="perc_perfect",
                                                                            plt_title="Percentage of Perfect Repeats by VRR Telomere Length",
                                                                            x_title="VRR Telomere Length (BP)", y_title="Percentage of Perfect Repeats",
                                                                            hover_tooltips=[("VRR Telomere Length", "@x"), ("Percentage Perfect Repeats", "@y")])
                            with Grid(columns=2):
                                EZChart(length_by_perfect, THEME)
                                EZChart(vrr_length_by_perfect, THEME)

                        with new_tabs.add_tab("Variant Repeats Distribution"):
                                # telo % GGTTAG {s<=1} hist
                                perfect_hist = report_utils.telo_length_hist(df["perc_variant"], binwidth=1, 
                                                                    binrange=[0, 100],
                                                                    plt_title="Telomeric Variant Repeat Composition",
                                                                    x_title="Perentage of Variant Repeats",
                                                                    y_title="Telomere Count")
                                # telo % GGTTAG boxplot
                                perfect_box = report_utils.create_boxplot(df=df,column_name="perc_variant", sample=sample,
                                                                            x_title = sample, y_title="Percentage of Telomere",
                                                                            plt_title="Telomeric Variant Repeat Composition")
                
                                # % GGTTAG by telomere quality
                                perc_by_quality = report_utils.scatterplot(data=df, x="telo_quality", y="perc_variant",
                                                                            plt_title="Telomere Quality by Composition",
                                                                            x_title="Telomere Quality", 
                                                                            y_title="Percentage of /Variant Repeats",
                                                                            hover_tooltips=[("Telomere Quality", "@x"), ("Percentage Variant Repeats", "@y")])
                                with Grid(columns=3):
                                    EZChart(perfect_hist, THEME)
                                    EZChart(perfect_box, THEME)
                                    EZChart(perc_by_quality, THEME)
                            
                        with new_tabs.add_tab("Variant Repeats by Telomere Length"):
                            # % GGTTAG by telomere length
                            length_by_perfect = report_utils.scatterplot(data=df, x="telo_length", y="perc_variant",
                                                                            plt_title="Percentage of Variant Repeats by Telomere Length",
                                                                            x_title="Telomere Length (BP)", y_title="Percentage of Variant Repeats",
                                                                            hover_tooltips=[("Telomere Length", "@x"), ("Percentage Variant Repeats", "@y")])
                            vrr_length_by_perfect = report_utils.scatterplot(data=df, x="vrr_telo_length", y="perc_variant",
                                                                            plt_title="Percentage of Variant Repeats by VRR Telomere Length",
                                                                            x_title="VRR Telomere Length (BP)", y_title="Percentage of Variant Repeats",
                                                                            hover_tooltips=[("VRR Telomere Length", "@x"), ("Percentage Variant Repeats", "@y")])
                            with Grid(columns=2):
                                EZChart(length_by_perfect, THEME)
                                EZChart(vrr_length_by_perfect, THEME)
                        with new_tabs.add_tab("Perfect vs Variant Repeats"):
                            scatter = report_utils.scatterplot(data=df, x="perc_perfect", y="perc_variant",
                                                                plt_title="Percentage of Telomere Composition",
                                                                x_title = "Percentage of Perfect Repeats",
                                                                y_title="Percentage of Variant Repeats")
                            boxplot = report_utils.create_multisample_boxplot(df=df, column_names=["perc_perfect", "perc_variant"],
                                                                                plt_title="Telomeric Composition Boxplot",
                                                                                y_title="Percentage of VRR Telomere", x_rotation=45,
                                                                                min_q=0, max_q=100)
                            # scatter plot showing these two?
                            # boxplot comparisons
                            with Grid(columns=2):
                                EZChart(scatter, THEME)
                                EZChart(boxplot, THEME)
            
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
    parser.add_argument("--sample_telo_stats", nargs="+", required=True)
    parser.add_argument("--run_telo_stats", required=True)
    parser.add_argument("--run_vrr_stats", required=True)
    parser.add_argument("--restriction_digest", required=True, nargs="+")
    parser.add_argument("--plot_vrr_length", required=True)
    parser.add_argument("--plot_telo_length", required=True)
    parser.add_argument("--strand_comparison", required=True)
    parser.add_argument("--detailed_stats", required=True)
    
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)

