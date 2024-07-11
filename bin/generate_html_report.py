#!/usr/bin/env python3

import sys
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import math
from jinja2 import Template
import argparse
import json


def main(args):

    params_file = open(args.params)
    params = params_file.read()
    params = json.loads(params)

    manifest_file = open(args.manifest)
    manifest = manifest_file.read()
    manifest = json.loads(manifest)

    versions = pd.read_table(args.versions, sep=",", header=None)
    versions.set_axis(["Software", "Version"], axis=1)

    versions_tables = go.Figure(data=[go.Table(
        header=dict(values=list(versions.columns),
                fill_color='paleturquoise',
                align='center'),
        cells=dict(values=versions.transpose().values.tolist(),
               fill_color='lavender',
               align='left'))
    ])

    jinja_data = {
        "pipeline_version":manifest["version"],
        "nextflow_version":manifest["nextflowVersion"], 
        "DOI":manifest["doi"],
        "homePage":manifest["homePage"],
        "author": manifest["author"],
        "run_name": params["run_name"],
        "input_file": params["input_file"],
        "out_dir": params["outdir"],
        "repeat": params["repeat"],
        "command": args.commandLine,
        "versions": versions_tables.to_html()
    }

    f = open(args.report, 'w')
    template_file = open(args.template_file, 'r')
    jinja_template = Template(template_file.read())
    f.write(jinja_template.render(jinja_data))
    f.close()
    template_file.close()   

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
    parser.add_argument("--sample_stats_retained", required=True) #works but only for simplex, not multiplex tested yet
    parser.add_argument("--sample_stats_filtered", required=True) #works but only for simplex, not multiplex tested yet
    
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)

