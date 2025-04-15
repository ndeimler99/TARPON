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
from dominate.tags import em, p, style
import ezcharts as ezc
from ezcharts.layout.base import IStyles, Snippet
from ezcharts.components.ezchart import EZChart
from ezcharts.components.fastcat import SeqSummary
#from ezcharts.components.reports.labs import LabsReport
from extra.report import LabsReport
#from report import BasicReport
from ezcharts.layout.snippets import DataTable
from ezcharts.layout.snippets import Grid
from ezcharts.layout.snippets import Tabs
import pandas as pd
#from extra.report_utils import BokehPlot as report_utils.BokehPlot
from ezcharts.components.reports.comp import ComponentReport
from bokeh.models import HoverTool, ColumnDataSource
import extra.report_utils as report_utils
import numpy as np
from bokeh.plotting import figure
import bokeh.palettes
from ezcharts.layout.util import css, render_template
from typing import Optional

THEME = 'epi2melabs'

def pass_fail_sample(row):

    if row["Number_of_Reads"] >= args.minimum_read_count:
        return "PASS"
    else:
        return "FAIL"
    
def get_nextflow_attributes(attribute_file):

    attribute_file = open(attribute_file)
    attribute = attribute_file.read()
    attribute = json.loads(attribute)
    return attribute

def main(args):

    args.strand_comparison = args.strand_comparison == "true"
    args.detailed_stats = args.detailed_stats == "true"
    args.minimum_read_count = int(args.minimum_read_count)

    params = get_nextflow_attributes(args.params)
    manifest = get_nextflow_attributes(args.manifest)

    versions = pd.read_table(args.versions, sep=",", header=None)
    versions.set_axis(["Software", "Version"], axis=1)

    sample_dict = {}
    for i in args.sample_stats_retained:
        sample_dict[i.split(".")[0]] = {}
        sample_dict[i.split(".")[0]]["retained_reads"] = i

    for i in args.sample_stats_removed:
        sample_dict[i.split(".")[0]]["removed_reads"] = i

    for i in args.sample_telo_stats:
        sample_dict[i.split(".")[0]]["telo_stats"] = i

    if args.restriction_digest[0] != "false":
        if i.endswith(".txt"):
            for i in args.restriction_digest:
                sample_dict[i.split(".")[0]]["digest"] = i

    vrr_summary_stats = pd.read_table(args.run_vrr_stats, sep="\t")

    report = LabsReport(
        f"Report for: {params['run_name']}", args.workflow_name,
        args.params, args.versions, manifest["version"], args.manifest)
    
    with report.add_section("Sequencing Stats", "Sequencing Stats"):
        p("""Statistics for the entire flow cell of non-demultiplexed data""")
        tabs = Tabs()
        with tabs.add_dropdown_menu("Retained Read Statistics", change_header=False):
            df = pd.read_table(args.run_stats_retained)
            df["file"] = df["file"].apply(lambda x: x.split("/")[-1].split(".bam")[0])
            if args.strand_comparison:
                g_strand = df[df.file.str.contains("g_strand")]
                c_strand = df[df.file.str.contains("c_strand")]
                g_strand["file"] = g_strand["file"].apply(lambda x: x.split(".g_strand")[0])
                c_strand["file"] = c_strand["file"].apply(lambda x: x.split(".c_strand")[0])
                i = df[(df.file.str.contains("strand"))].index
                df = df.drop(i)
            with tabs.add_dropdown_tab('Read Count'):
                plt = report_utils.barplot(data=df, x="file", y="num_seqs", 
                                           x_title="Pipeline Step", x_rotation=45,
                                           y_title="Number of Retained_Reads",
                                           plt_title="Number of Reads Retained at Each Step")
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Number of Reads Retained", "@num_seqs")]
                EZChart(plt,THEME) 

                if args.strand_comparison:
                    g_plt = report_utils.barplot(data=g_strand, x="file", y="num_seqs", 
                                           x_title="Pipeline Step", x_rotation=45,
                                           y_title="Number of Retained_Reads",
                                           plt_title="G Strand")
                    hover = g_plt._fig.select(dict(type=HoverTool))
                    hover.tooltips = [("Number of Reads Retained", "@num_seqs")]

                    c_plt = report_utils.barplot(data=c_strand, x="file", y="num_seqs", 
                                           x_title="Pipeline Step", x_rotation=45,
                                           y_title="Number of Retained_Reads",
                                           plt_title="C Strand")
                    hover = c_plt._fig.select(dict(type=HoverTool))
                    hover.tooltips = [("Number of Reads Retained", "@num_seqs")]
                    with Grid(columns=2):
                        EZChart(g_plt, THEME)
                        EZChart(c_plt, THEME)

            with tabs.add_dropdown_tab('Stats Table'):
                DataTable.from_pandas(df, use_index=False)
                if args.strand_comparison:
                    with Grid(columns=2):
                        DataTable.from_pandas(g_strand, use_index=False)
                        DataTable.from_pandas(c_strand, use_index=False)

            with tabs.add_dropdown_tab('Read Length'):
                plt = report_utils.seqkit_stats_boxplot_length(df, x="file", 
                                                               y_title="Read Length (BP)", x_title="Pipeline Step", 
                                                               x_rotation=45,
                                                               plt_title="Read Length of Reads Retained at Each Step")
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Length", "@mean_length"),
                                  ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                   ("Min Length", "@min_length"), ("Max Length", "@max_length")]
                EZChart(plt, THEME)

                if args.strand_comparison:
                    g_plt = report_utils.seqkit_stats_boxplot_length(g_strand, x="file", 
                                                               y_title="Read Length (BP)", x_title="Pipeline Step", 
                                                               x_rotation=45,
                                                               plt_title="G Strand")
                    hover = g_plt._fig.select(dict(type=HoverTool))
                    hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Length", "@mean_length"),
                                  ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                   ("Min Length", "@min_length"), ("Max Length", "@max_length")]
                    c_plt = report_utils.seqkit_stats_boxplot_length(data=c_strand, x="file", 
                                                               y_title="Read Length (BP)", x_title="Pipeline Step", 
                                                               x_rotation=45,
                                                               plt_title="C Strand")
                    hover = c_plt._fig.select(dict(type=HoverTool))
                    hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Length", "@mean_length"),
                                  ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                   ("Min Length", "@min_length"), ("Max Length", "@max_length")]
                    with Grid(columns=2):
                        EZChart(g_plt, THEME)
                        EZChart(c_plt, THEME)

            with tabs.add_dropdown_tab('Read Quality'):
                df = pd.read_table(args.run_stats_retained)
                df["file"] = df["file"].apply(lambda x: x.split("/")[-1].split(".bam")[0])
                if args.strand_comparison:
                    g_strand = df[df['file'].str.contains("g_strand")]
                    c_strand = df[df['file'].str.contains("c_strand")]
                    g_strand["file"] = g_strand["file"].apply(lambda x: x.split(".g_strand")[0])
                    c_strand["file"] = c_strand["file"].apply(lambda x: x.split(".c_strand")[0])
                    i = df[(df.file.str.contains("strand"))].index
                    df = df.drop(i)
                plt = report_utils.quality_boxplot_from_quantiles(data=df, x="file",
                                           x_title="Pipeline Step", y_title="Quality Score",
                                           plt_title="Quality of Reads Retained at Each Step",
                                           x_rotation=45)
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Quality", "@mean_quality"),
                                  ("Q1", "@Q1_qual"),("Q2", "@Q2_qual"),("Q3", "@Q3_qual"),
                                   ("Min Quality", "@min_quality"), ("Max Quality", "@max_quality")]
                EZChart(plt,THEME)

                if args.strand_comparison:
                    g_plt = report_utils.quality_boxplot_from_quantiles(data=g_strand, x="file",
                                            x_title="Pipeline Step", y_title="Quality Score",
                                            plt_title="G Strand",
                                            x_rotation=45)
                    hover = g_plt._fig.select(dict(type=HoverTool))
                    hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Quality", "@mean_quality"),
                                  ("Q1", "@Q1_qual"),("Q2", "@Q2_qual"),("Q3", "@Q3_qual"),
                                   ("Min Quality", "@min_quality"), ("Max Quality", "@max_quality")]
                    c_plt = report_utils.quality_boxplot_from_quantiles(data=c_strand, x="file",
                                            x_title="Pipeline Step", y_title="Quality Score",
                                            plt_title="C Strand",
                                            x_rotation=45)
                    hover = c_plt._fig.select(dict(type=HoverTool))
                    hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Quality", "@mean_quality"),
                                  ("Q1", "@Q1_qual"),("Q2", "@Q2_qual"),("Q3", "@Q3_qual"),
                                   ("Min Quality", "@min_quality"), ("Max Quality", "@max_quality")]
                    with Grid(columns=2):
                        EZChart(g_plt, THEME)
                        EZChart(c_plt, THEME)

        with tabs.add_dropdown_menu("Removed Read Statistics", change_header=False):
            df = pd.read_table(args.run_stats_removed)
            df["file"] = df["file"].apply(lambda x: x.split("/")[-1].split(".bam")[0])
            
            if args.strand_comparison:
                g_strand = df[df['file'].str.contains("g_strand")]
                c_strand = df[df['file'].str.contains("c_strand")]
                g_strand["file"] = g_strand["file"].apply(lambda x: x.split(".g_strand")[0])
                c_strand["file"] = c_strand["file"].apply(lambda x: x.split(".c_strand")[0])
                i = df[(df.file.str.contains("strand"))].index
                df = df.drop(i)
            with tabs.add_dropdown_tab('Read Count'):
                plt=report_utils.barplot(data=df, x="file", y="num_seqs",
                                plt_title="Number of Reads Removed at Each Step",
                                x_title="Pipeline Step", x_rotation=45,
                                y_title="Number of Removed Reads")
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Number of Reads Removed", "@num_seqs")]
                # what is this function doing
                EZChart(plt,THEME)

                if args.strand_comparison:
                    g_plt = report_utils.barplot(data=g_strand, x="file", y="num_seqs", 
                                           x_title="Pipeline Step", x_rotation=45,
                                           y_title="Number of Removed Reads",
                                           plt_title="G Strand")
                    hover = g_plt._fig.select(dict(type=HoverTool))
                    hover.tooltips = [("Number of Reads Removed", "@num_seqs")]

                    c_plt = report_utils.barplot(data=c_strand, x="file", y="num_seqs", 
                                           x_title="Pipeline Step", x_rotation=45,
                                           y_title="Number of Removed Reads",
                                           plt_title="C Strand")
                    hover = c_plt._fig.select(dict(type=HoverTool))
                    hover.tooltips = [("Number of Reads Removed", "@num_seqs")]
                    with Grid(columns=2):
                        EZChart(g_plt, THEME)
                        EZChart(c_plt, THEME)

            with tabs.add_dropdown_tab('Stats Table'):
                DataTable.from_pandas(df, use_index=False)
                if args.strand_comparison:
                    with Grid(columns=2):
                        DataTable.from_pandas(g_strand, use_index=False)
                        DataTable.from_pandas(c_strand, use_index=False)

            with tabs.add_dropdown_tab('Read Length'):
                plt = report_utils.seqkit_stats_boxplot_length(df, x="file",
                                                               y_title="Read Length (BP)", 
                                                               x_title="Pipeline Step", x_rotation=45,
                                                               plt_title="Read Length of Reads Removed at Each Step")
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Length", "@mean_length"),
                                  ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                   ("Min Length", "@min_length"), ("Max Length", "@max_length")]
                # what is this function doing?
                EZChart(plt, THEME)

                if args.strand_comparison:
                    g_plt = report_utils.seqkit_stats_boxplot_length(g_strand, x="file", 
                                                               y_title="Read Length (BP)", x_title="Pipeline Step", 
                                                               x_rotation=45,
                                                               plt_title="G Strand")
                    hover = g_plt._fig.select(dict(type=HoverTool))
                    hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Length", "@mean_length"),
                                  ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                   ("Min Length", "@min_length"), ("Max Length", "@max_length")]
                    c_plt = report_utils.seqkit_stats_boxplot_length(data=c_strand, x="file", 
                                                               y_title="Read Length (BP)", x_title="Pipeline Step", 
                                                               x_rotation=45,
                                                               plt_title="C Strand")
                    hover = c_plt._fig.select(dict(type=HoverTool))
                    hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Length", "@mean_length"),
                                  ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                   ("Min Length", "@min_length"), ("Max Length", "@max_length")]
                    with Grid(columns=2):
                        EZChart(g_plt, THEME)
                        EZChart(c_plt, THEME)


            with tabs.add_dropdown_tab('Read Quality'):
                df = pd.read_table(args.run_stats_removed)
                df["file"] = df["file"].apply(lambda x: x.split("/")[-1].split(".bam")[0])

                if args.strand_comparison:
                    g_strand = df[df['file'].str.contains("g_strand")]
                    c_strand = df[df['file'].str.contains("c_strand")]
                    g_strand["file"] = g_strand["file"].apply(lambda x: x.split(".g_strand")[0])
                    c_strand["file"] = c_strand["file"].apply(lambda x: x.split(".c_strand")[0])
                    i = df[(df.file.str.contains("strand"))].index
                    df = df.drop(i)

                plt = report_utils.quality_boxplot_from_quantiles(data=df, x="file",
                                           x_title="Pipeline Step", y_title="Quality Score",
                                           plt_title="Quality of Reads Removed at Each Step",
                                           x_rotation=45)
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Quality", "@mean_quality"),
                                  ("Q1", "@Q1_qual"),("Q2", "@Q2_qual"),("Q3", "@Q3_qual"),
                                   ("Min Length", "@min_quality"), ("Max Length", "@max_quality")]
                EZChart(plt,THEME)

                if args.strand_comparison:
                    g_plt = report_utils.quality_boxplot_from_quantiles(data=g_strand, x="file",
                                            x_title="Pipeline Step", y_title="Quality Score",
                                            plt_title="G Strand",
                                            x_rotation=45)
                    hover = g_plt._fig.select(dict(type=HoverTool))
                    hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Quality", "@mean_quality"),
                                  ("Q1", "@Q1_qual"),("Q2", "@Q2_qual"),("Q3", "@Q3_qual"),
                                   ("Min Length", "@min_quality"), ("Max Length", "@max_quality")]
                    c_plt = report_utils.quality_boxplot_from_quantiles(data=c_strand, x="file",
                                            x_title="Pipeline Step", y_title="Quality Score",
                                            plt_title="C Strand",
                                            x_rotation=45)
                    hover = c_plt._fig.select(dict(type=HoverTool))
                    hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Quality", "@mean_quality"),
                                  ("Q1", "@Q1_qual"),("Q2", "@Q2_qual"),("Q3", "@Q3_qual"),
                                   ("Min Length", "@min_quality"), ("Max Length", "@max_quality")]
                    with Grid(columns=2):
                        EZChart(g_plt, THEME)
                        EZChart(c_plt, THEME)
        
    with report.add_section("Sample Comparison", "Sample Comparison"):
        p("""Between Sample Comparisons of Telomere Length""")
        tabs = Tabs()
        with tabs.add_tab("Number of Telomeric Reads"):
            #barplot using numseqs
       
            df = vrr_summary_stats

            df = df.set_index("Sample_ID").loc[sorted(sample_dict.keys())].reset_index()
            plt=report_utils.barplot(data=df, x="Sample_ID", y="Number_of_Reads",
                                     plt_title="Number of Telomeric Reads per Sample",
                                     x_title="Sample", x_rotation=45,
                                     y_title="Number of Telomeres")
            hover = plt._fig.select(dict(type=HoverTool))

            hover.tooltips = [("Sample", "@Sample_ID"), ("Number of Reads", "@Number_of_Reads"),("Mean_Telomere_Length", "@Mean_VRR_Telomere_Length")]

            EZChart(plt, THEME)
            ####

        with tabs.add_tab("VRR Stats"):
            df = vrr_summary_stats
            df = df.set_index("Sample_ID").loc[sorted(sample_dict.keys())].reset_index()
            df["STATUS"] = df.apply(pass_fail_sample, axis=1)
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
                                
                    if args.strand_comparison:
                        with new_tabs.add_tab("VRR Telo Length Strand Comparison"):
                            with Grid(columns=3):
                                with Grid(columns=1):
                                    plt = report_utils.telo_length_hist_by_strand([df[df["strand"] == "G"]["vrr_telo_length"], df[df["strand"] == "C"]["vrr_telo_length"]], labels=["G strand", "C strand"],
                                                                                binwidth=200, binrange=[0,max(df["vrr_telo_length"])+200],
                                                                                plt_title="Telomere Length G-Strand",
                                                                                x_title="Telomere Length (BP)", y_title="Percentage of Sequences")
                                    EZChart(plt, THEME)

                                telo_hist = report_utils.create_boxplot_by_strand(df, "vrr_telo_length", plt_title="Telo Length by Strand", x_title="Strand", y_title="Telomere Length (bp)")
                                EZChart(telo_hist, THEME)
                                    

                                bins = [i*1000 for i in range(0,11)]
                                bins.append(100000)
                                telo_bar_df_g = np.histogram(df[df["strand"]=="G"]["vrr_telo_length"], bins=bins)
                                telo_bar_df_g = pd.DataFrame(list(zip(telo_bar_df_g[1], telo_bar_df_g[0])), columns=["bin_start", "bin_size"])
                                telo_bar_df_g["bin_start"] = telo_bar_df_g["bin_start"].astype("string")
                                telo_bar_df_g["sample"] = "G"
                                telo_bar_df_g["bin_size"] = telo_bar_df_g["bin_size"] / sum(telo_bar_df_g["bin_size"]) * 100

                                telo_bar_df_c = np.histogram(df[df["strand"]=="C"]["vrr_telo_length"], bins=bins)
                                telo_bar_df_c = pd.DataFrame(list(zip(telo_bar_df_c[1], telo_bar_df_c[0])), columns=["bin_start", "bin_size"])
                                telo_bar_df_c["bin_start"] = telo_bar_df_c["bin_start"].astype("string")
                                telo_bar_df_c["sample"] = "C"
                                telo_bar_df_c["bin_size"] = telo_bar_df_c["bin_size"] / sum(telo_bar_df_c["bin_size"]) * 100


                                telo_bar_plot = report_utils.telo_barplot(data=pd.concat([telo_bar_df_g, telo_bar_df_c]), 
                                                        x="sample", x_rotation=45, x_title="Strand",
                                                        y="bin_size", y_title="Percentage of Reads",
                                                        hue="bin_start", dodge=False,
                                                        order=["G", "C"], 
                                                        palette=bokeh.palettes.Category20[11],
                                                        legend_loc="right", legend_orientation="vertical",
                                                        plt_title="Telomere Length by Strand Binned Bar Plot")
                                EZChart(telo_bar_plot, THEME)

                                
                with tabs.add_dropdown_tab("{} Retained Reads".format(sample)):
                    new_tabs = Tabs()
                    df = pd.read_table(sample_dict[sample]["retained_reads"], sep="\t")
                    df["file"] = df["file"].apply(lambda x: x.split("/")[-1].split(".bam")[0])
                    if args.strand_comparison:
                        g_strand = df[df.file.str.contains("g_strand")]
                        c_strand = df[df.file.str.contains("c_strand")]
                        g_strand["file"] = g_strand["file"].apply(lambda x: x.split(".g_strand")[0])
                        c_strand["file"] = c_strand["file"].apply(lambda x: x.split(".c_strand")[0])
                        i = df[(df.file.str.contains("strand"))].index
                        df = df.drop(i)

                    with new_tabs.add_tab("Read Count"):
                        #bar plot
                        plt = report_utils.barplot(data=df, x="file", y="num_seqs",
                                                   plt_title="Number of Reads Retained at Each Step",
                                                   x_title="Pipeline Step", x_rotation=45, y_title="Number of Retained Reads")
                        hover = plt._fig.select(dict(type=HoverTool))
                        hover.tooltips = [("Number of Reads Retained", "@num_seqs")]
                        EZChart(plt,THEME)

                        if args.strand_comparison:
                            g_plt = report_utils.barplot(data=g_strand, x="file", y="num_seqs", 
                                                x_title="Pipeline Step", x_rotation=45,
                                                y_title="Number of Retained_Reads",
                                                plt_title="G Strand")
                            hover = g_plt._fig.select(dict(type=HoverTool))
                            hover.tooltips = [("Number of Reads Retained", "@num_seqs")]

                            c_plt = report_utils.barplot(data=c_strand, x="file", y="num_seqs", 
                                                x_title="Pipeline Step", x_rotation=45,
                                                y_title="Number of Retained_Reads",
                                                plt_title="C Strand")
                            hover = c_plt._fig.select(dict(type=HoverTool))
                            hover.tooltips = [("Number of Reads Retained", "@num_seqs")]
                            with Grid(columns=2):
                                EZChart(g_plt, THEME)
                                EZChart(c_plt, THEME)

                    #stats
                    with new_tabs.add_tab("Read Stats"):
                        DataTable.from_pandas(df, use_index=False)
                        if args.strand_comparison:
                            with Grid(columns=2):
                                DataTable.from_pandas(g_strand, use_index=False)
                                DataTable.from_pandas(c_strand, use_index=False)
                    #length boxplot
                    with new_tabs.add_tab("Read Length"):
                        plt = report_utils.seqkit_stats_boxplot_length(df, x="file",
                                                                       plt_title="Read Length of Retained Reads",
                                                                       x_title="Pipeline Step", x_rotation=45,
                                                                       y_title="Read Length (BP)")
                        hover = plt._fig.select(dict(type=HoverTool))
                        hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                        ("Avg Length", "@mean_length"),
                                        ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                        ("Min Length", "@min_length"), ("Max Length", "@max_length")]
                        EZChart(plt, THEME)

                        if args.strand_comparison:
                            g_plt = report_utils.seqkit_stats_boxplot_length(g_strand, x="file", 
                                                                    y_title="Read Length (BP)", x_title="Pipeline Step", 
                                                                    x_rotation=45,
                                                                    plt_title="G Strand")
                            hover = g_plt._fig.select(dict(type=HoverTool))
                            hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                        ("Avg Length", "@mean_length"),
                                        ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                        ("Min Length", "@min_lengtg"), ("Max Length", "@max_length")]
                            c_plt = report_utils.seqkit_stats_boxplot_length(data=c_strand, x="file", 
                                                                    y_title="Read Length (BP)", x_title="Pipeline Step", 
                                                                    x_rotation=45,
                                                                    plt_title="C Strand")
                            hover = c_plt._fig.select(dict(type=HoverTool))
                            hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                        ("Avg Length", "@mean_length"),
                                        ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                        ("Min Length", "@min_lengtg"), ("Max Length", "@max_lengtg")]
                            with Grid(columns=2):
                                EZChart(g_plt, THEME)
                                EZChart(c_plt, THEME)

                    #quality plot
                    with new_tabs.add_tab("Read Quality"):
                        df = pd.read_table(sample_dict[sample]["retained_reads"])
                        df["file"] = df["file"].apply(lambda x: x.split("/")[-1].split(".bam")[0])
                        if args.strand_comparison:
                            g_strand = df[df.file.str.contains("g_strand")]
                            c_strand = df[df.file.str.contains("c_strand")]
                            g_strand["file"] = g_strand["file"].apply(lambda x: x.split(".g_strand")[0])
                            c_strand["file"] = c_strand["file"].apply(lambda x: x.split(".c_strand")[0])
                            df = df.drop(i)
                        plt = report_utils.quality_boxplot_from_quantiles(data=df, x="file",
                                           x_title="Pipeline Step", y_title="Quality Score",
                                           plt_title="Quality of Reads Retained at Each Step",
                                           x_rotation=45)
                        hover = plt._fig.select(dict(type=HoverTool))
                        hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Quality", "@mean_quality"),
                                  ("Q1", "@Q1_qual"),("Q2", "@Q2_qual"),("Q3", "@Q3_qual"),
                                   ("Min Quality", "@min_quality"), ("Max Quality", "@max_quality")]
                        EZChart(plt,THEME)

                        if args.strand_comparison:

                            g_plt = report_utils.quality_boxplot_from_quantiles(data=g_strand, x="file",
                                                    x_title="Pipeline Step", y_title="Quality Score",
                                                    plt_title="G Strand",
                                                    x_rotation=45)
                            hover = g_plt._fig.select(dict(type=HoverTool))
                            hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Quality", "@mean_quality"),
                                  ("Q1", "@Q1_qual"),("Q2", "@Q2_qual"),("Q3", "@Q3_qual"),
                                   ("Min Quality", "@min_quality"), ("Max Quality", "@max_quality")]
                            c_plt = report_utils.quality_boxplot_from_quantiles(data=c_strand, x="file",
                                                    x_title="Pipeline Step", y_title="Quality Score",
                                                    plt_title="C Strand",
                                                    x_rotation=45)
                            hover = c_plt._fig.select(dict(type=HoverTool))
                            hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Quality", "@mean_quality"),
                                  ("Q1", "@Q1_qual"),("Q2", "@Q2_qual"),("Q3", "@Q3_qual"),
                                   ("Min Quality", "@min_quality"), ("Max Quality", "@max_quality")]
                            with Grid(columns=2):
                                EZChart(g_plt, THEME)
                                EZChart(c_plt, THEME)

                with tabs.add_dropdown_tab("{} Removed Reads".format(sample)):
                    new_tabs = Tabs()
                    df = pd.read_table(sample_dict[sample]["removed_reads"], sep="\t")
                    df["file"] = df["file"].apply(lambda x: x.split("/")[-1].split(".bam")[0])
                    if args.strand_comparison:
                        g_strand = df[df.file.str.contains("g_strand")]
                        c_strand = df[df.file.str.contains("c_strand")]
                        g_strand["file"] = g_strand["file"].apply(lambda x: x.split(".g_strand")[0])
                        c_strand["file"] = c_strand["file"].apply(lambda x: x.split(".c_strand")[0])
                        i = df[(df.file.str.contains("strand"))].index
                        df = df.drop(i)
                    with new_tabs.add_tab("Read Count"):
                        #bar plot
                        plt = report_utils.barplot(data=df, x="file", y="num_seqs", 
                                                   plt_title="Number of Reads Removed at Each Step",
                                                   x_title="Pipeline Step", x_rotation=45,
                                                   y_title="Number of Retained Reads")
                        hover = plt._fig.select(dict(type=HoverTool))
                        hover.tooltips = [("Number of Reads Removed", "@num_seqs")]
                        EZChart(plt,THEME)
                        if args.strand_comparison:
                            g_plt = report_utils.barplot(data=g_strand, x="file", y="num_seqs", 
                                                x_title="Pipeline Step", x_rotation=45,
                                                y_title="Number of Removed Reads",
                                                plt_title="G Strand")
                            hover = g_plt._fig.select(dict(type=HoverTool))
                            hover.tooltips = [("Number of Reads Removed", "@num_seqs")]

                            c_plt = report_utils.barplot(data=c_strand, x="file", y="num_seqs", 
                                                x_title="Pipeline Step", x_rotation=45,
                                                y_title="Number of Removed Reads",
                                                plt_title="C Strand")
                            hover = c_plt._fig.select(dict(type=HoverTool))
                            hover.tooltips = [("Number of Reads Removed", "@num_seqs")]
                            with Grid(columns=2):
                                EZChart(g_plt, THEME)
                                EZChart(c_plt, THEME)

                    #stats
                    with new_tabs.add_tab("Read Stats"):
                        DataTable.from_pandas(df, use_index=False)
                        if args.strand_comparison:
                            with Grid(columns=2):
                                DataTable.from_pandas(g_strand, use_index=False)
                                DataTable.from_pandas(c_strand, use_index=False)
                    #length boxplot
                    with new_tabs.add_tab("Read Length"):
                        plt = report_utils.seqkit_stats_boxplot_length(df, x="file",
                                                                       y_title="Read Length (BP)",
                                                                       x_title="Pipeline Step",
                                                                       x_rotation=45)
                        hover = plt._fig.select(dict(type=HoverTool))
                        hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                        ("Avg Length", "@mean_length"),
                                        ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                        ("Min Length", "@min_lengtg"), ("Max Length", "@max_length")]
                        EZChart(plt, THEME)
                        if args.strand_comparison:
                            g_plt = report_utils.seqkit_stats_boxplot_length(g_strand, x="file", 
                                                                    y_title="Read Length (BP)", x_title="Pipeline Step", 
                                                                    x_rotation=45,
                                                                    plt_title="G Strand")
                            hover = g_plt._fig.select(dict(type=HoverTool))
                            hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                        ("Avg Length", "@mean_length"),
                                        ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                        ("Min Length", "@min_length"), ("Max Length", "@max_length")]
                            c_plt = report_utils.seqkit_stats_boxplot_length(data=c_strand, x="file", 
                                                                    y_title="Read Length (BP)", x_title="Pipeline Step", 
                                                                    x_rotation=45,
                                                                    plt_title="C Strand")
                            hover = c_plt._fig.select(dict(type=HoverTool))
                            hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                        ("Avg Length", "@mean_length"),
                                        ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                        ("Min Length", "@min_length"), ("Max Length", "@max_length")]
                            with Grid(columns=2):
                                EZChart(g_plt, THEME)
                                EZChart(c_plt, THEME)
                    #quality plot
                    with new_tabs.add_tab("Read Quality"):
                        df = pd.read_table(sample_dict[sample]["removed_reads"])
                        df["file"] = df["file"].apply(lambda x: x.split("/")[-1].split(".bam")[0])
                        if args.strand_comparison:
                            g_strand = df[df.file.str.contains("g_strand")]
                            c_strand = df[df.file.str.contains("c_strand")]
                            g_strand["file"] = g_strand["file"].apply(lambda x: x.split(".g_strand")[0])
                            c_strand["file"] = c_strand["file"].apply(lambda x: x.split(".c_strand")[0])
                            i = df[(df.file.str.contains("strand"))].index
                            df = df.drop(i)
                        plt = report_utils.quality_boxplot_from_quantiles(data=df, x="file",
                                           x_title="Pipeline Step", y_title="Quality Score",
                                           plt_title="Quality of Reads Retained at Each Step",
                                           x_rotation=45)
                        hover = plt._fig.select(dict(type=HoverTool))
                        hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Quality", "@mean_quality"),
                                  ("Q1", "@Q1_qual"),("Q2", "@Q2_qual"),("Q3", "@Q3_qual"),
                                   ("Min Length", "@min_quality"), ("Max Length", "@max_quality")]
                        EZChart(plt,THEME)
                        if args.strand_comparison:
            
                            g_plt = report_utils.quality_boxplot_from_quantiles(data=g_strand, x="file",
                                                    x_title="Pipeline Step", y_title="Quality Score",
                                                    plt_title="G Strand",
                                                    x_rotation=45)
                            hover = g_plt._fig.select(dict(type=HoverTool))
                            hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Quality", "@mean_quality"),
                                  ("Q1", "@Q1_qual"),("Q2", "@Q2_qual"),("Q3", "@Q3_qual"),
                                   ("Min Length", "@min_quality"), ("Max Length", "@max_quality")]
                            c_plt = report_utils.quality_boxplot_from_quantiles(data=c_strand, x="file",
                                                    x_title="Pipeline Step", y_title="Quality Score",
                                                    plt_title="C Strand",
                                                    x_rotation=45)
                            hover = c_plt._fig.select(dict(type=HoverTool))
                            hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Quality", "@mean_quality"),
                                  ("Q1", "@Q1_qual"),("Q2", "@Q2_qual"),("Q3", "@Q3_qual"),
                                   ("Min Length", "@min_quality"), ("Max Length", "@max_quality")]
                            with Grid(columns=2):
                                EZChart(g_plt, THEME)
                                EZChart(c_plt, THEME)
                with tabs.add_dropdown_tab("{} Telomere Composition".format(sample)):
                    new_tabs = Tabs()
                    stats_df = pd.read_table(sample_dict[sample]["telo_stats"], sep="\t")
                    with new_tabs.add_tab("Histogram and Boxplot"):
                        if args.mutant == "false":
                            histogram_repeat_freq = report_utils.repeat_freq_histogram(stats_df[["wt_composition","one_nucl_variant_composition"]],
                                                                                    binwidth=1, binrange=[0,100],
                                                                                    x_title="% of Telomere", y_title="Number of Telomeric Reads",
                                                                                    plt_title="Repeat Frequency")
                            
                            
                    
                            # boxplot of repeat frequencies
                            repeat_freq_boxplot = report_utils.create_multisample_boxplot(stats_df, ["wt_composition", "one_nucl_variant_composition"], 
                                                                                        0, 100, plt_title="Repeat Frequency",
                                                                                        y_title = "Percentage of Telomere")
                            # telomere length barplot w/ colors
                            with Grid(columns=2):
                                EZChart(histogram_repeat_freq, THEME)
                                EZChart(repeat_freq_boxplot, THEME)
                        else:
                            # histogram of repeat frequencies
                            histogram_repeat_freq = report_utils.repeat_freq_histogram(stats_df[["wt_composition","mutant_composition", "one_nucl_variant_composition"]],
                                                                                    binwidth=1, binrange=[0,100],
                                                                                    x_title="% of Telomere", y_title="Number of Telomeric Reads",
                                                                                    plt_title="Repeat Frequency")
                            
                            
                    
                            # boxplot of repeat frequencies
                            repeat_freq_boxplot = report_utils.create_multisample_boxplot(stats_df, ["wt_composition", "mutant_composition", "one_nucl_variant_composition"], 
                                                                                        0, 100, plt_title="Repeat Frequency",
                                                                                        y_title = "Percentage of Telomere")
                            # telomere length barplot w/ colors
                            with Grid(columns=2):
                                EZChart(histogram_repeat_freq, THEME)
                                EZChart(repeat_freq_boxplot, THEME)
                        # telomere length barplot w/ colors
                    with new_tabs.add_tab("Telomere Length by Composition Bar Plot"):
                        stats_df = stats_df.sort_values("vrr_telo_length")
                        if args.mutant == "false":
                            plot = report_utils.colored_telo_length_barplot(data=stats_df, x_title="VRR Telomere Length (bp)", y_title="Read ID", repeat=args.repeat)
                        else:
                            plot = report_utils.colored_telo_length_barplot(data=stats_df, mutant=args.mutant, x_title="VRR Telomere Length (bp)", y_title="Read ID", repeat=args.repeat)
                        EZChart(plot, THEME)
                if args.restriction_digest[0] != "false":
                        with tabs.add_dropdown_tab("{} Restriction Digest".format(sample)):
                            stats_df = pd.read_table(sample_dict[sample]["digest"], sep="\t")
                            stats_df["file"] = stats_df["file"].apply(lambda x: x.split(".reheaded.bam")[0])
                            DataTable.from_pandas(stats_df, use_index=False)
    
    if args.detailed_stats:
        with report.add_section("Detailed Analysis", "Detailed Analysis"):
            tabs = Tabs()
            for sample in sorted(list(sample_dict.keys())):
                df = pd.read_table(sample_dict[sample]["telo_stats"], sep="\t")
                with tabs.add_dropdown_menu(sample, change_header=False):
                    with tabs.add_dropdown_tab("Sequencing Quality"):
                        new_tabs = Tabs()
                        with new_tabs.add_tab("Read vs. Telo Quality Comparison"):
                            # read quality hist
                            read_qual = report_utils.telo_length_hist(df["read_qual"], binwidth=1, 
                                                                    binrange=[0, max(df["read_qual"] + 1)],
                                                                    plt_title="Read Quality",
                                                                    x_title="Quality Score",
                                                                    y_title="Read Count")
            
                            # telo quality hist
                            telo_qual = report_utils.telo_length_hist(df["telo_qual"], binwidth=1, 
                                                                    binrange=[0, max(df["telo_qual"] + 1)],
                                                                    plt_title="Telomere Quality",
                                                                    x_title="Quality Score",
                                                                    y_title="Telomere Count")
                            # read quality vs telo scatterplot
                            read_telo_scatter = report_utils.scatterplot(data=df, x="telo_qual", y="read_qual",
                                                                        plt_title="Telomere vs Read Sequencing Quality",
                                                                        x_title= "Telomere Quality", y_title="Read Quality",
                                                                       hover_tooltips=[("Telomere Quality", "@x"), ("Read Quality", "@y")])
                            
                            
                            boxplot = report_utils.create_multisample_boxplot(df=df, column_names=["telo_qual", "read_qual"], 
                                                                                min_q = 15, max_q = 50, plt_title="Seq. Quality",
                                                                                y_title="Quality Score", x_rotation=45,
                                                                                x_labels={"telo_qual":"Telomere Quality", "read_qual":"Read Quality"})
                            
                            with Grid(columns=4):
                                EZChart(read_qual, THEME)
                                EZChart(telo_qual, THEME)            
                                EZChart(boxplot, THEME)
                                EZChart(read_telo_scatter, THEME)
                                
                        with new_tabs.add_tab("Quality vs VRR Length"):
                                # read quality by vrr length
                            read_by_telo = report_utils.scatterplot(data=df, x="read_qual", y="vrr_telo_length",
                                                                        plt_title="Read Quality by VRR Telo Length",
                                                                        x_title="Read Quality", y_title="VRR Telomere Length (BP)",
                                                                        hover_tooltips=[("Read Quality", "@x"), ("VRR Telomere Length", "@y")])
                                # telo quality by telo length
                            telo_by_telo = report_utils.scatterplot(data=df, x="telo_qual", y="vrr_telo_length",
                                                                        plt_title="Telomere Quality by VRR Telo Length",
                                                                        x_title="Telomere Quality", y_title="VRR Telomere Length (BP)",
                                                                        hover_tooltips=[("Telomere Quality", "@x"), ("VRR Telomere Length", "@y")])
                                
                            with Grid(columns=2):
                                EZChart(read_by_telo, THEME)
                                EZChart(telo_by_telo, THEME)
                        if args.strand_comparison:
                            with new_tabs.add_tab("Strand Comparison"):
                                
                                read_quality = report_utils.create_boxplot_by_strand(df, "read_qual", 
                                    plt_title="Read Quality by Strand", x_title="Strand", y_title="Read Quality")
                                telo_quality = report_utils.create_boxplot_by_strand(df, "telo_qual", 
                                    plt_title="Telo Quality by Strand", x_title="Strand", y_title="Telomere Quality")

                                with Grid(columns=2):
                                    EZChart(read_quality, THEME)
                                    EZChart(telo_quality, THEME)
                                # read quality by strand
                                # telo quality by strand

                                pass
                    with tabs.add_dropdown_tab("Telomeric Composition"):
                        new_tabs = Tabs()
                        with new_tabs.add_tab("Perfect Repeats Distribution"):
                                # telo % GGTTAG hist
                                perfect_hist = report_utils.telo_length_hist(df["wt_composition"], binwidth=1, 
                                                                    binrange=[0, 100],
                                                                    plt_title="Telomeric Perfect Repeat Composition",
                                                                    x_title="Perentage of Perfect Repeats",
                                                                    y_title="Telomere Count")
                                # telo % GGTTAG boxplot
                                perfect_box = report_utils.create_boxplot(df=df,column_name="wt_composition", sample=sample,
                                                                            x_title = sample, y_title="Percentage of Telomere",
                                                                            plt_title="Telomeric Perfect Repeat Composition")
                
                                # % GGTTAG by telomere quality
                                perc_by_quality = report_utils.scatterplot(data=df, x="telo_qual", y="wt_composition",
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
                            length_by_perfect = report_utils.scatterplot(data=df, x="telo_length", y="wt_composition",
                                                                            plt_title="Percentage of Perfect Repeats by Telomere Length",
                                                                            x_title="Telomere Length (BP)", y_title="Percentage of Perfect Repeats",
                                                                            hover_tooltips=[("Telomere Length", "@x"), ("Percentage Perfect Repeats", "@y")])
                            vrr_length_by_perfect = report_utils.scatterplot(data=df, x="vrr_telo_length", y="wt_composition",
                                                                            plt_title="Percentage of Perfect Repeats by VRR Telomere Length",
                                                                            x_title="VRR Telomere Length (BP)", y_title="Percentage of Perfect Repeats",
                                                                            hover_tooltips=[("VRR Telomere Length", "@x"), ("Percentage Perfect Repeats", "@y")])
                            with Grid(columns=2):
                                EZChart(length_by_perfect, THEME)
                                EZChart(vrr_length_by_perfect, THEME)

                        with new_tabs.add_tab("Variant Repeats Distribution"):
                                # telo % GGTTAG {s<=1} hist
                                perfect_hist = report_utils.telo_length_hist(df["one_nucl_variant_composition"], binwidth=1, 
                                                                    binrange=[0, 100],
                                                                    plt_title="Telomeric Variant Repeat Composition",
                                                                    x_title="Perentage of Variant Repeats",
                                                                    y_title="Telomere Count")
                                # telo % GGTTAG boxplot
                                perfect_box = report_utils.create_boxplot(df=df,column_name="one_nucl_variant_composition", sample=sample,
                                                                            x_title = sample, y_title="Percentage of Telomere",
                                                                            plt_title="Telomeric Variant Repeat Composition")
                
                                # % GGTTAG by telomere quality
                                perc_by_quality = report_utils.scatterplot(data=df, x="telo_qual", y="one_nucl_variant_composition",
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
                            length_by_perfect = report_utils.scatterplot(data=df, x="telo_length", y="one_nucl_variant_composition",
                                                                            plt_title="Percentage of Variant Repeats by Telomere Length",
                                                                            x_title="Telomere Length (BP)", y_title="Percentage of Variant Repeats",
                                                                            hover_tooltips=[("Telomere Length", "@x"), ("Percentage Variant Repeats", "@y")])
                            vrr_length_by_perfect = report_utils.scatterplot(data=df, x="vrr_telo_length", y="one_nucl_variant_composition",
                                                                            plt_title="Percentage of Variant Repeats by VRR Telomere Length",
                                                                            x_title="VRR Telomere Length (BP)", y_title="Percentage of Variant Repeats",
                                                                            hover_tooltips=[("VRR Telomere Length", "@x"), ("Percentage Variant Repeats", "@y")])
                            with Grid(columns=2):
                                EZChart(length_by_perfect, THEME)
                                EZChart(vrr_length_by_perfect, THEME)
                        with new_tabs.add_tab("Perfect vs Variant Repeats"):
                            scatter = report_utils.scatterplot(data=df, x="wt_composition", y="one_nucl_variant_composition",
                                                                plt_title="Percentage of Telomere Composition",
                                                                x_title = "Percentage of Perfect Repeats",
                                                                y_title="Percentage of Variant Repeats")
                            boxplot = report_utils.create_multisample_boxplot(df=df, column_names=["wt_composition", "one_nucl_variant_composition"],
                                                                                plt_title="Telomeric Composition Boxplot",
                                                                                y_title="Percentage of VRR Telomere", x_rotation=45,
                                                                                min_q=0, max_q=100)
                            # scatter plot showing these two?
                            # boxplot comparisons
                            with Grid(columns=2):
                                EZChart(scatter, THEME)
                                EZChart(boxplot, THEME)
                        if args.strand_comparison:
                            with new_tabs.add_tab("Strand Comparison"):
                                # perfect repeat percentage by strand
                                perf_repeat = report_utils.create_boxplot_by_strand(df, "wt_composition", 
                                                                                    plt_title="Perfect Repeat Composition", x_title="Strand", 
                                                                                    y_title="Percentage Perfect Repeats")

                                # imperfect repeat percentage by strand
                                imperf_repeat = report_utils.create_boxplot_by_strand(df, "one_nucl_variant_composition", 
                                                                                    plt_title="Telomere-Like Repeat Composition", x_title="Strand", 
                                                                                    y_title="Percentage Telomere-Like Repeats")

                                with Grid(columns=2):
                                    EZChart(perf_repeat, THEME)
                                    EZChart(imperf_repeat, THEME)
    
    if args.mutant != "false":
        with report.add_section("Mutant Repeat Analysis", "Mutant Repeat Analysis"):
            pass

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
    parser.add_argument("--minimum_read_count", required=True)
    parser.add_argument("--run_stats_retained", required=True) #works but need to figure out how to get combined sample stats here
    parser.add_argument("--run_stats_removed", required=True) #works but need to figure out how to get combined sample stats here
    parser.add_argument("--sample_stats_retained", nargs='+', required=True) #works but only for simplex, not multiplex tested yet
    parser.add_argument("--sample_stats_removed", nargs='+', required=True) #works but only for simplex, not multiplex tested yet
    parser.add_argument("--sample_telo_stats", nargs="+", required=True)
    parser.add_argument("--run_vrr_stats", required=True)
    parser.add_argument("--restriction_digest", required=True, nargs="+")
    parser.add_argument("--strand_comparison", required=True)
    parser.add_argument("--detailed_stats", required=True)
    parser.add_argument("--mutant", required=True)
    parser.add_argument("--mutant_analysis_repeat_distribution", nargs="+", required=True)
    parser.add_argument("--mutant_analysis_processivity", nargs="+", required=True)
    parser.add_argument("--repeat", required=True)
    
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)

