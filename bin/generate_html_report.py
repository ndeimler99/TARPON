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


#  Sample Comparison Telomere Length Barchart – add hover
#  Sample Comparison Telomere Length Boxplot – add hover
#        a.	Use custom boxplot function that takes in dataframe containing precalculated Quartile values – similar to seqkit boxplot function

#  Individual Sample – histogram add hover
#  Individual Sample boxplot add hover
#  Individual Sample scatterplot add better hover
#  Individual Sample histogram and boxplot add red line denoting mean
#  Individual sample boxplot and barchart remove x axis tick labels
#  Add Manifest tablel similar to params table
#  Report page – remove provided by Oxford Nanopore Technologies
#  Report page – remove About footer – or at least modify text
#  Report page – remove Epi2Me logo and replace with TArPON cartoon



THEME = 'epi2melabs'

def get_nextflow_attributes(attribute_file):

    attribute_file = open(attribute_file)
    attribute = attribute_file.read()
    attribute = json.loads(attribute)
    return attribute

def main(args):

    params = get_nextflow_attributes(args.params)
    # params_file = open(args.params)
    # params = params_file.read()
    # params = json.loads(params)

    manifest = get_nextflow_attributes(args.manifest)
    # manifest_file = open(args.manifest)
    # manifest = manifest_file.read()
    # manifest = json.loads(manifest)


    versions = pd.read_table(args.versions, sep=",", header=None)
    versions.set_axis(["Software", "Version"], axis=1)

    sample_dict = {}
    for i in args.sample_stats_retained:
        sample_dict[i.split(".")[0]] = {}
        sample_dict[i.split(".")[0]]["retained_reads"] = i

    for i in args.sample_stats_filtered:
        sample_dict[i.split(".")[0]]["filtered_reads"] = i

    run_telo_stats = pd.DataFrame()
    for i in args.sample_telo_stats:
        sample_dict[i.split(".")[0]]["telo_stats"] = i
        df = pd.read_table(i)
        df["sample"] = i.split(".")[0]
        run_telo_stats = pd.concat([run_telo_stats, df], ignore_index=True, sort=False)

    if args.restriction_digest[0] != "false":
        if i.endswith(".txt"):
            for i in args.restriction_digest:
                sample_dict[i.split(".")[0]]["digest"] = i
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
            with tabs.add_dropdown_tab('Read Count'):
                #convert to custom barplot function
                plt = report_utils.barplot(data=df, x="file", y="num_seqs", 
                                  x_title="Pipeline Step", x_rotation=45,
                                  y_title="Number of Retained_Reads",
                                  plt_title="Number of reads Retained at Each Step")
                
                #plt._fig.title.text = "Number of Reads Retained at Each Step"
                #plt._fig.xaxis.axis_label = "Pipeline Step"
                #plt._fig.yaxis.axis_label = "Number of Retained Reads"
                #plt._fig.xaxis.major_label_orientation = 45
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Number of Reads Retained", "@top")]
                # what is this function doing?
                EZChart(plt,THEME)
            with tabs.add_dropdown_tab('Stats Table'):
                # this is fine
                DataTable.from_pandas(df, use_index=False)
            with tabs.add_dropdown_tab('Read Length'):
                # add yaxis and axis label to function call
                plt = report_utils.seqkit_stats_boxplot_length(df, x="file")
                plt._fig.yaxis.axis_label = "Read Length (BP)"
                plt._fig.xaxis.axis_label = "Pipeline Step"
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Length", "@avg_len"),
                                  ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                   ("Min Length", "@min_len"), ("Max Length", "@max_len")]
                EZChart(plt, THEME)

            with tabs.add_dropdown_tab('Read Quality'):
                # add title, axis, etc to custom function
                plt = report_utils.barplot(data=df, x="file", y="AvgQual")
                plt._fig.title.text = "Average Quality of Reads Retained at Each Step"
                plt._fig.xaxis.axis_label = "Pipeline Step"
                plt._fig.yaxis.axis_label = "Average Quality"
                plt._fig.xaxis.major_label_orientation = 45
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("File", "@file"), ("Average Quality", "@AvgQual"), ("Read Count", "@num_seqs")]
                EZChart(plt,THEME)

        with tabs.add_dropdown_menu("Filtered Read Statistics", change_header=False):
            df = pd.read_table(args.run_stats_filtered)
            df["file"] = df["file"].apply(lambda x: x.split("/")[-1].split(".fastq")[0])
            df = df.drop(columns=["N50.1", "format", "sum_len", "sum_gap"])

            with tabs.add_dropdown_tab('Read Count'):
                # convert to custom function with axis calls and title in function
                plt=ezc.barplot(data=df, x="file", y="num_seqs")
                plt._fig.title.text = "Number of Reads Retained at Each Step"
                plt._fig.xaxis.axis_label = "Pipeline Step"
                plt._fig.yaxis.axis_label = "Number of Filtered Reads"
                plt._fig.xaxis.major_label_orientation = 45
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Number of Reads Filtered", "@top")]
                # what is this function doing
                EZChart(plt,THEME)
            with tabs.add_dropdown_tab('Stats Table'):
                # this is fine
                DataTable.from_pandas(df, use_index=False)
            with tabs.add_dropdown_tab('Read Length'):
                # add title, axis, etc to custom function
                plt = report_utils.seqkit_stats_boxplot_length(df, x="file")
                plt._fig.yaxis.axis_label = "Read Length (BP)"
                plt._fig.xaxis.axis_label = "Pipeline Step"
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                  ("Avg Length", "@avg_len"),
                                  ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                   ("Min Length", "@min_len"), ("Max Length", "@max_len")]
                # what is this function doing?
                EZChart(plt, THEME)

            with tabs.add_dropdown_tab('Read Quality'):
                # add xaxis, yaxis, title to this function
                plt=report_utils.barplot(data=df, x="file", y="AvgQual")
                plt._fig.title.text = "Average Quality of Reads Filtered at Each Step"
                plt._fig.xaxis.axis_label = "Pipeline Step"
                plt._fig.yaxis.axis_label = "Average Quality"
                plt._fig.xaxis.major_label_orientation = 45
                hover = plt._fig.select(dict(type=HoverTool))
                hover.tooltips = [("File", "@file"), ("Average Quality", "@AvgQual"), ("Read Count", "@num_seqs")]
                # what is this function doing
                EZChart(plt,THEME)

    with report.add_section("Sample Comparison", "Sample Comparison"):
        p("""Between Sample Comparisons of Telomere Length""")
        tabs = Tabs()
        with tabs.add_tab("Number of Telomeric Reads"):
            #barplot using numseqs
            df = pd.read_table(args.run_telo_stats_table, sep="\t")
            df = df.set_index("Sample_ID").loc[sorted(sample_dict.keys())].reset_index()
            palette = bokeh.palettes.Category20[20]
            df["color"] = palette[0:len(df["Sample_ID"])]
            #### New Test
            source = ColumnDataSource(df)
            plt = BokehPlot()
            
            plt._fig = figure(x_range=[str(x) for x in df["Sample_ID"]], y_range=(0, max(df["Number_of_Reads"])*1.1),
                              tools="hover,crosshair,pan,box_zoom,zoom_in,zoom_out,reset,save")
            plt._fig.vbar(x="Sample_ID", top="Number_of_Reads", width=0.9, source=source, color="color")
            plt._fig.toolbar.logo = None
            axis_color = "#677884"
            axis_width = 1.2
            axis_font = "normal"
            plt._fig.yaxis.axis_line_color = axis_color
            plt._fig.yaxis.axis_line_width = axis_width
            plt._fig.yaxis.axis_label_text_font_style = axis_font
            plt._fig.xaxis.axis_line_color = axis_color
            plt._fig.xaxis.axis_line_width = axis_width
            plt._fig.xaxis.axis_label_text_font_style = axis_font
            plt._fig.yaxis.minor_tick_line_color = None
            plt._fig.xaxis.major_tick_line_color = axis_color
            plt._fig.yaxis.major_tick_line_color = axis_color
            plt._fig.title.text_font_size = "18px"
            plt._fig.xaxis.axis_label = "Sample"
            plt._fig.yaxis.axis_label = "Number of Telomeres"
            hover = plt._fig.select(dict(type=HoverTool))
            hover.tooltips = [("Sample", "@Sample_ID"), ("Number of Reads", "@Number_of_Reads"),("Mean_Telomere_Length", "@Mean_Telomere_Length")]
            EZChart(plt, THEME)
            ####

            #plt = ezc.barplot(data=df, x="Sample_ID", y="Number of Reads", order=sorted(list(sample_dict.keys())))
            #plt._fig.xaxis.axis_label = "Sample ID"
            #hover = plt._fig.select(dict(type=HoverTool))
            #hover.tooltips = [("Sample", "@Sample_ID"), ("Read Count", "@top")]
            #EZChart(plt, THEME)
        with tabs.add_tab("Telo Stats"):
            df = pd.read_table(args.run_telo_stats_table, sep="\t")
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
            

            telo_bar_plot = ezc.barplot(data=master_df, x="sample", y="bin_size", hue="bin_start", dodge=False, order=sorted(list(sample_dict.keys())), palette=bokeh.palettes.Category20[11])
            
            telo_bar_plot._fig.yaxis.axis_label = "Percentage of Telomeres"
            telo_bar_plot._fig.xaxis.axis_label = "Sample"
            hover = telo_bar_plot._fig.select(dict(type=HoverTool))
            hover.tooltips = [("Sample", "@sample"),("% of Telomeres", "@bin_Size")]
            EZChart(telo_bar_plot, THEME)
        with tabs.add_tab("Telomere Length Boxplot"):
            df = pd.read_table(args.run_telo_stats, sep=",")
            plt = ezc.boxplot(data=df, x="sample", y="telo_length", order=sorted(list(sample_dict.keys())))
            hover = plt._fig.select(dict(type=HoverTool))
            hover.tooltips = [("Sample", "@Sample_ID"),("Number of Reads", "@Number_of_Reads")]
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
                    with new_tabs.add_tab("Telo Length Analysis (n={})".format(len(df['telo_length']))):
                        #telo_length_box = ezc.boxplot()
                        # telo_small_bins = np.histogram(df["telo_length"], bins=[i*100 for i in range(0, round(max(df['telo_length'])/100) + 1)])
                        # telo_small_bins = pd.DataFrame(list(zip(telo_small_bins[1], telo_small_bins[0])), columns=["bin_start", "bin_size"])
                        # telo_small_bins["bin_min"] = telo_small_bins["bin_start"] + 1
                        # telo_small_bins["bin_max"] = telo_small_bins["bin_start"] + 100
                        # #print(telo_small_bins)
                        # telo_length_hist = report_utils.telo_hist_plot(data=telo_small_bins, x="bin_start", y="bin_size", color="blue")
                        # hover = telo_length_hist._fig.select(dict(type=HoverTool))
                        # hover.tooltips = [("Bin Size", "@bin_min - @bin_max"), ("Count", "@top")]
                        # telo_length_hist._fig.yaxis.axis_label = "Telomere Count"
                        # telo_length_hist._fig.xaxis.axis_label = "Telomere Length (BP)"
                        telo_length_hist = ezc.histplot(df["telo_length"], binwidth=100, binrange=[0,max(df["telo_length"])+100])
                        telo_length_hist.xAxis.name = "Telo Length (BP)"
                        telo_length_hist.yAxis.name = "Read Count"
                        telo_length_boxplot = report_utils.create_boxplot(df=df, column_name="telo_length")
                        telo_length_boxplot._fig.yaxis.axis_label = "Telomere Length"
                        telo_length_boxplot._fig.xaxis.axis_label = sample
                        # telo_length_boxplot._fig.title.text = "{} Telomere Length Distribution".format(sample)
                        bins = [i*1000 for i in range(0,11)]
                        bins.append(100000)
                        telo_bar_df = np.histogram(df["telo_length"], bins=bins)
                        telo_bar_df = pd.DataFrame(list(zip(telo_bar_df[1], telo_bar_df[0])), columns=["bin_start", "bin_size"])
                        telo_bar_df["bin_start"] = telo_bar_df["bin_start"].astype("string")
                        telo_bar_df["sample"] = sample
                        telo_bar_df["bin_size"] = telo_bar_df["bin_size"] / sum(telo_bar_df["bin_size"])
                        #telo_bar_df["bin_label"] = "{} - {} bp".format(telo_bar_df["bin_start"] + 1, telo_bar_df["bin_start"]+1000)
                        telo_bar_plot = ezc.barplot(data=telo_bar_df, x="sample", y="bin_size", hue="bin_start", dodge=False, palette=bokeh.palettes.Category20[11])
                        telo_bar_plot._fig.yaxis.axis_label = "Percentage of Telomeres"
                        telo_bar_plot._fig.xaxis.axis_label = sample


                        #telo_length_bar = "tmp"
                        with Grid(columns=3):
                            EZChart(telo_length_hist, THEME)
                            EZChart(telo_length_boxplot, THEME)
                            EZChart(telo_bar_plot, THEME)
                    with new_tabs.add_tab("Telo vs Read Length"):
                        telo_length_hist = ezc.histplot(df["telo_length"], binwidth=100, binrange=[0,max(df["telo_length"])+100])
                        telo_length_hist.xAxis.name = "Telo Length (BP)"
                        telo_length_hist.yAxis.name = "Read Count"
                        read_length_hist = ezc.histplot(df["read_len"], binwidth=100, binrange=[0, max(df["read_len"])+100])
                        read_length_hist.xAxis.name = "Read Length (BP)"
                        read_length_hist.yAxis.name = "Read Count"
                        byscatter = ezc.scatterplot(data=df, x="read_len", y="telo_length")
                        byscatter.xAxis.name = "Read Length"
                        byscatter.yAxis.name = "Telomere Length"
                        with Grid(columns=3):
                            #telo length hist, read length_hist, scatter plot
                            EZChart(telo_length_hist, THEME)
                            EZChart(read_length_hist, THEME)
                            EZChart(byscatter, THEME)
                
                with tabs.add_dropdown_tab("{} Retained Reads".format(sample)):
                    new_tabs = Tabs()
                    df = pd.read_table(sample_dict[sample]["retained_reads"], sep="\t")
                    df["file"] = df["file"].apply(lambda x: x.split("/")[-1].split(".fastq")[0])
                    df = df.drop(columns=["N50.1", "format", "sum_len", "sum_gap"])
                    with new_tabs.add_tab("Read Length"):
                        #bar plot
                        plt = ezc.barplot(data=df, x="file", y="num_seqs")
                        plt._fig.title.text = "Number of Reads Retained at Each Step"
                        plt._fig.xaxis.axis_label = "Pipeline Step"
                        plt._fig.yaxis.axis_label = "Number of Retained Reads"
                        plt._fig.xaxis.major_label_orientation = 45
                        hover = plt._fig.select(dict(type=HoverTool))
                        hover.tooltips = [("Number of Reads Retained", "@top")]
                        EZChart(plt,THEME)
                    #stats
                    with new_tabs.add_tab("Read Stats"):
                        DataTable.from_pandas(df, use_index=False)
                    #length boxplot
                    with new_tabs.add_tab("Read Length"):
                        plt = report_utils.seqkit_stats_boxplot_length(df, x="file")
                        plt._fig.yaxis.axis_label = "Read Length (BP)"
                        plt._fig.xaxis.axis_label = "Pipeline Step"
                        hover = plt._fig.select(dict(type=HoverTool))
                        hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                        ("Avg Length", "@avg_len"),
                                        ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                        ("Min Length", "@min_len"), ("Max Length", "@max_len")]
                        EZChart(plt, THEME)
                    #quality plot
                    with new_tabs.add_tab("Read Quality"):
                        plt=report_utils.barplot(data=df, x="file", y="AvgQual")
                        plt._fig.title.text = "Average Quality of Reads Retained at Each Step"
                        plt._fig.xaxis.axis_label = "Pipeline Step"
                        plt._fig.yaxis.axis_label = "Average Quality"
                        plt._fig.xaxis.major_label_orientation = 45
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
                        plt = ezc.barplot(data=df, x="file", y="num_seqs")
                        plt._fig.title.text = "Number of Reads Filtered at Each Step"
                        plt._fig.xaxis.axis_label = "Pipeline Step"
                        plt._fig.yaxis.axis_label = "Number of Retained Reads"
                        plt._fig.xaxis.major_label_orientation = 45
                        hover = plt._fig.select(dict(type=HoverTool))
                        hover.tooltips = [("Number of Reads Filtered", "@top")]
                        EZChart(plt,THEME)
                    #stats
                    with new_tabs.add_tab("Read Stats"):
                        DataTable.from_pandas(df, use_index=False)
                    #length boxplot
                    with new_tabs.add_tab("Read Length"):
                        plt = report_utils.seqkit_stats_boxplot_length(df, x="file")
                        plt._fig.yaxis.axis_label = "Read Length (BP)"
                        plt._fig.xaxis.axis_label = "Pipeline Step"
                        hover = plt._fig.select(dict(type=HoverTool))
                        hover.tooltips = [("Sample", "@file"), ("Read Count", "@num_seqs"),
                                        ("Avg Length", "@avg_len"),
                                        ("Q1", "@Q1"),("Q2", "@Q2"),("Q3", "@Q3"),
                                        ("Min Length", "@min_len"), ("Max Length", "@max_len")]
                        EZChart(plt, THEME)
                    #quality plot
                    with new_tabs.add_tab("Read Quality"):
                        plt=report_utils.barplot(data=df, x="file", y="AvgQual")
                        plt._fig.title.text = "Average Quality of Reads Filtered at Each Step"
                        plt._fig.xaxis.axis_label = "Pipeline Step"
                        plt._fig.yaxis.axis_label = "Average Quality"
                        plt._fig.xaxis.major_label_orientation = 45
                        hover = plt._fig.select(dict(type=HoverTool))
                        hover.tooltips = [("File", "@file"), ("Average Quality", "@AvgQual"), ("Read Count", "@num_seqs")]
                        EZChart(plt,THEME)
                
                if args.restriction_digest[0] != "false":
                    with tabs.add_dropdown_tab("{} Restriction Digest".format(sample)):
                        stats_df = pd.read_table(sample_dict[sample]["digest"], sep="\t")
                        stats_df = stats_df.drop(columns=["N50.1", "format", "sum_len", "sum_gap"])
                        DataTable.from_pandas(stats_df, use_index=False)
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
    parser.add_argument("--sample_telo_stats", nargs="+", required=True)
    parser.add_argument("--run_telo_stats_table", required=True)
    parser.add_argument("--run_telo_stats", required=True)
    parser.add_argument("--restriction_digest", default=False, nargs="+")
    #    parser.add_argument("--restriction_digest_analysis", required=False, default="test")
    
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)

