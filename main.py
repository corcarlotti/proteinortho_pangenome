# to-do:
# remove all multiple occurrences of same functions and rename
# 2 functions for bar and pie plot

# merge all scripts to a command line executable programm

import pandas as pd
import random
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import glob

from add_desc_to_id import add_desc
from protoeinortho_processing import *
from proteinortho_visualization import *
from gene_contribution import *

#text="This is a tool for analysing the proteinortho output in regard of a pangenome analysis."

parser = argparse.ArgumentParser()#description=text)
parser.version = "1.0"
parser.add_argument("-i", "--inputfile", nargs="+", help="path to tabular proteinortho output file(s) \nrequired for options -A/-G/-N/-V")#, action="store_true")
parser.add_argument("-d", "--description", help="path to tabular proteinortho description file")#, action="store_false")
parser.add_argument("-N", "--NAMES", help="add gene names to IDs for the tabular input file", action="store_true")
parser.add_argument("-A", "--ANALYSE", help="tabular pangenome analysis", action="store_true")
parser.add_argument("-V", "--VISUALIZE", help="visualization of the pangenome composition as pie or bar plot", action="store_true")
parser.add_argument("-G", "--GENE_CONTRIBUTION", help="box plot of gene contribution with each genome", action="store_true")
parser.add_argument("-v", "--version", help="show program version", action="version")
#args = parser.parse_args() #read arguments from command line
args, unknown = parser.parse_known_args()


# Check for function paramaters

# input file
def check_inputfile():
    if args.inputfile:
        return(True)
    else:
        print("Please enter an input file with the parameters -i/--inputfile. Check the usage with the -h/--help option.")
        return(False)

# description file
def check_description():
    if args.description:
        return(True)
    else:
        print("Please enter a description file with the parameters -d/--description. Check the usage with the -h/--help option.")
        return(False)

# NAMES
def names():
    for i in args.inputfile:
        outputfile='.'.join(i.split(".")[:-1]+["named.tsv"])
        add_desc(i, args.description, outputfile)

# ANALYSE
def analyse():
    result_df = create_df()
    for tab in args.inputfile:
        table = pd.read_csv(tab, sep='\t', lineterminator='\n')
        desc = pd.read_csv(args.description, sep='\t', lineterminator='\n', header=None, names=["IDs", "Description"])
        filled_df = fill_df(table, desc, result_df, tab)
    outputfile='.'.join(tab.split(".")[:1]+["analysis.tsv"])
    filled_df.to_csv(outputfile)

# VISUALIZE
def visualize():
    # bar/pie?
    for tab in args.inputfile:
        table = pd.read_csv(tab, sep='\t', lineterminator='\n')
        desc = pd.read_csv(args.description, sep='\t', lineterminator='\n', header=None, names=["IDs", "Description"])
        plotfile = '.'.join(tab.split(".")[:-1] + ["plots.png"])
        pangenome_visualize(table, desc, plotfile)

# GENE CONTRIBUTION
def gene_contr():
    for tab in args.inputfile:
        table = pd.read_csv(tab, sep='\t', lineterminator='\n')
        table_description = pd.read_csv(args.description, sep='\t', lineterminator='\n', header=None, names=["IDs", "Description"])

        # build base of strain and identifier
        base = build_base(table)

        # calculate unique genes
        result_dict = get_gene_count(table_description)
        unique_dict = calculate_unique_genes(table, base, result_dict)

        # calculate and visualize data
        rep_number = [100, 1000]
        out_name = '.'.join(tab.split(".")[:-1] + ["boxplot.png"])
        visualize_box_plots(rep_number, table, table_description, base, unique_dict, out_name)

def check_function_arguments():
    if args.NAMES:
        names()
    elif args.ANALYSE:
        analyse()
    elif args.VISUALIZE:
        visualize()
    elif args.GENE_CONTRIBUTION:
        gene_contr()
    else:
        print("Please enter function argument. Check the usage with the -h/--help option.")
        exit()

def main():

    i_file = check_inputfile()
    d_file = check_description()
    if i_file == True & d_file == True:
        check_function_arguments()

if __name__ == "__main__":
    main()
