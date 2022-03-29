#!/usr/bin/env python
# coding: utf-8

import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from reportlab.lib import colors
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from dna_features_viewer import load_record 
from dna_features_viewer import BiopythonTranslator
from dna_features_viewer import BiopythonTranslator

print(sys.argv[0]) # prints python_script.py
print("input path for circuitgff file:\t", sys.argv[1]) # prints var1
print("input path for depth file of sense fragments:\t",sys.argv[2])
print("input path for depth file of antisense fragments:\t",sys.argv[3])
print("output path for coverage depth plot with circuit annotation:\t",sys.argv[4])
print("sample name for plot shown as:\t", sys.argv[5])
print("ploting requested figure...\n")

def plotdepwithanno(gfffile, sendepth, antisendepth, plotfile, samp):
    fig, (ax1, ax2, ax3) = plt.subplots(
 3, 1, figsize=(30,6), sharex=True,gridspec_kw={"height_ratios": [6,1, 1]}
)
    class CustomTranslator(BiopythonTranslator):
 # Label fields indicates the order in which annotations fields are
 # considered to determine the feature's label
    #label_fields = ["label", "note", "name", "gene"]
        label_fields = ["label", "Name"]
        def compute_feature_legend_text(self, feature):
            return feature.type
    

        def compute_feature_color(self, feature):
            return {
            "gene": "yellow",
             "promoter": "#ffd383", # light orange
             "ribozyme": "red",
             "terminator": "#fbf3f6", # pink
             "transcript": "#d1e9f1", # light blue
             "promoter_unit": "darkblue",
             }[feature.type]


        def compute_feature_box_color(self, feature):
            return "white"

        def compute_feature_box_linewidth(self, feature):
            return 0
    
    translator = CustomTranslator()
    graphic_record = translator.translate_record(gfffile)
    graphic_record.plot(ax=ax1, with_ruler=False, strand_in_label_threshold=7)
    graphic_record.plot_legend(ax=ax1, loc=1, frameon=False)
    #ax.figure.savefig("A_linear_plot.svg", bbox_inches="tight") 
    def plotdepth(file):
        depth = pd.read_table(file, sep="\t", header = None)
        depth = depth.rename(columns={0: "chr", 1: "loc", 2: "dep"})
        circuit = depth[depth['chr'] == "0x58v50"]
    
        return circuit["dep"]
    
    def plotloc(file):
        depth = pd.read_table(file, sep="\t", header = None)
        depth = depth.rename(columns={0: "chr", 1: "loc", 2: "dep"})
        circuit = depth[depth['chr'] == "0x58v50"]
    
        return circuit["loc"]

    ax1.set_title("circuit_annotations", fontsize=10)
    ax2.plot(plotloc(sendepth), plotdepth(sendepth))
    t2 = "_".join([samp, "sense", "coverage","depth","by","site"])
    ax2.set_title(t2, fontsize=10)
    x3.plot(plotloc(antisendepth), plotdepth(antisendepth),'tab:green')
    t3 = "_".join([samp, "antisense", "coverage","depth","by","site"])
    ax3.set_title(t3, fontsize=10)
    
    return fig.savefig(plotfile, bbox_inches="tight")


plotdepwithanno(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

print("job completed, please check it out in", sys.argv[4])

