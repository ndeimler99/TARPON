import sys
import os
import pandas as pd
import plotly.express as px
import math
from jinja2 import Template

output_dir = sys.argv[1]
pipeline_version = sys.argv[2]
nextflow_version = sys.argv[3]

telo_stats = pd.read_table("{}/FINAL_TELO/telomeric_stats.txt".format(output_dir))
bin_width=100
telo_bin_number = int(math.ceil(telo_stats['telo_length'].max() / 100.0)) * 2
telo_histogram = px.histogram(telo_stats, x="telo_length", nbins=telo_bin_number, labels={'telo_length': 'Telomere Length (Base Pairs)'})
telo_histogram.update_layout(yaxis_title="Number of Telomeric Reads",  margin=dict(l=6, r=6, t=6, b=6)) 
telo_histogram.update_xaxes(title_standoff=1)


telo_lengths_for_binning = [500, 1500, 2500, 3500, 4500, 5500, 6500, 7500, 8500, 9500, 10500]
telo_binned = [ min(telo_lengths_for_binning, key=lambda x:abs(x-length)) for length in telo_stats['telo_length']]

telo_counts = []
for x in telo_lengths_for_binning:
    telo_counts.append(telo_binned.count(x))

sample = ["sample" for i in telo_counts]
telo_counts = [x / sum(telo_counts) * 100 for x in telo_counts]

telo_lengths_for_binning = ["1-1,000", "1,001-2,000", "2,001-3,000", "3,001-4,000", "4,001-5,000", "5,001-6,000", "6,001-7,000", "7,001-8,000", "8,001-9,000", "9,001-10,000", "10,001+"]
telo_bar_df = pd.DataFrame({"lengths":telo_lengths_for_binning, "counts":telo_counts, "sample":sample})

telo_bar = px.bar(telo_bar_df, x="sample", y="counts", color="lengths", 
                  labels={"counts":"Percentage of Telomeric Sequences (%)", "sample":"", "lengths":"Telomere Length (BP)"},
                  color_discrete_map={
                    "1-1,000":"#FF68A1", 
                    "1,001-2,000":"#FF61CC", 
                    "2,001-3,000":"#C77CFF", 
                    "3,001-4,000":"#8494FF", 
                    "4,001-5,000":"#00A9FF", 
                    "5,001-6,000":"#00BFC4", 
                    "6,001-7,000":"#00BE67", 
                    "7,001-8,000":"#0CB702", 
                    "8,001-9,000":"#ABA300", 
                    "9,001-10,000":"#E68613", 
                    "10,001+":"#F8766D"
                    }
                )

telo_bar.update_xaxes(showticklabels=False)
telo_bar.update_layout(margin=dict(l=6, r=6, t=6, b=2), legend_y=0.5) 


read_counts_table = pd.read_table("{}/STATS/telomeric_stats.txt".format(output_dir))
read_counts_table = read_counts_table.filter(items=['file', 'num_seqs', 'avg_len', 'AvgQual'])

read_counts_table['file'] = read_counts_table['file'].apply(lambda x: x.split("/")[-1].split(".fastq")[0])
read_counts_table = read_counts_table.loc[read_counts_table['file'].isin(["input", "telomeric"])]


read_counts = read_counts_table.to_html(index=False).replace('<table border="1" class="dataframe">','<table border="1" class="table table-striped">')



retained_reads = pd.read_table("{}/STATS/telomeric_stats.txt".format(output_dir))
retained_reads['file'] = retained_reads['file'].apply(lambda x: x.split("/")[-1].split(".fastq")[0])
#retained_reads = retained_reads.loc[retained_reads['file'] != 'input']
retained_reads['file'] = pd.Categorical(retained_reads['file'], ["input", "putative_reads", "putative_reads.c_g_filtered", "subtelo", "adaptor", "telomeric"])
retained_reads = retained_reads.sort_values("file")

retained_reads_plot = px.bar(retained_reads, x="file", y='num_seqs', 
                             range_y=[0, sorted(list(set(retained_reads['num_seqs'])))[-2] + 1000], text_auto=True,
                             labels={'num_seqs':"Number of Reads", 'file':"Pipeline Process Output"})
retained_reads_plot.update_layout(margin=dict(l=6, r=6, t=6, b=6)) 


filtered_reads = pd.read_table("{}/STATS/filtered_stats.txt".format(output_dir))
filtered_reads['file'] = filtered_reads['file'].apply(lambda x: x.split("/")[-1].split(".fastq")[0])
filtered_reads['file'] = pd.Categorical(filtered_reads['file'], ["input", "non_telomeric", "20_80_removed_reads", "subtelo_filtered", "adaptor_filtered", 
                                                                 "no_telomere_start", "below_telo_%_threshold"])
filtered_reads = filtered_reads.sort_values("file")

filtered_reads_plot = px.bar(filtered_reads, x="file", y='num_seqs', 
                             range_y=[0, sorted(list(set(filtered_reads['num_seqs'])))[-3] + 1000], text_auto=True,
                             labels={'num_seqs':"Number of Reads", 'file':"Pipeline Process Output"})
filtered_reads_plot.update_layout(margin=dict(l=6, r=6, t=6, b=6)) 

if sys.argv[6]:
    mode = "C Strand Only"
else:
    mode = "Duplex Approach"

jinja_data = {
            "telo_histogram_plot":telo_histogram.to_html(full_html=False),
            "telo_bar_plot":telo_bar.to_html(full_html=False),
            "n":len(telo_stats['telo_length']),
            "pipeline_version":pipeline_version,
            "nextflow_version":nextflow_version,
            "read_counts":read_counts,
            "retained_reads_plot":retained_reads_plot.to_html(full_html=False),
            "filtered_reads_plot":filtered_reads_plot.to_html(full_html=False),
            "input_file":sys.argv[4],
            "out_dir":sys.argv[6],
            "repeat": sys.argv[5],
            "mode": mode,
            }

f = open(sys.argv[8], 'w')
template_file = open(sys.argv[7], 'r')
jinja_template = Template(template_file.read())
f.write(jinja_template.render(jinja_data))
f.close()
template_file.close()
