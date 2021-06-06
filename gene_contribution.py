# python script for calculation of gene contribution with adding genomes to the pangenome in a random order of set iterations.
# boxplot is created and a set of two iterations are compared in two subplots

import pandas as pd
import random
import numpy as np
import matplotlib.pyplot as plt
import sys, getopt


#build base with dict of all included strain names and their identifier
def build_base(tab):
    strain = []    # key and value list initialization
    identifier = []
    for i in tab.columns[3:]: # for every species name
        strain.append(i)      # save species name 
        ident = pd.Series.to_string(tab[tab['# Species'] == (tab.shape[1] - 3)][i][:1], index=False)    # select row with id
        identifier.append(ident.split('_')[0][0:])  # save cut id after '_'
        base_dict = dict(zip(strain, identifier))   # build a dict 
    return (base_dict)

# count all total (unique, core, acc) genes of one genome
def get_gene_count(table_desc):
    # group by cut ids in description and create dict
    table_desc["IDs"] = table_desc["IDs"].str.rsplit('_').str[0] 
    allgenes_dict = dict(table_desc.groupby(["IDs"]).size()) 
    return (allgenes_dict)

# count cluster included genes per genome 
def genes_in_column(table, strain, identifier):
    #count identifier in table column
    counter = 0                         
    for row in table[strain]:               
        counter += row.count(identifier)
    return (counter)

# fast version of counting unique genes (not clustered in table)
def calculate_unique_genes(table, value_dict, counter_dict):
    ug_id, ug_count = [], []    # initialize
    for key, val in value_dict.items(): # for every id-strain pair from dict
        # for every genome: total genes - used genes 
        used_genes = genes_in_column(table, key, val)
        total_genes = counter_dict.get(val) 
        ug_id.append(val) # save genome id
        not_used = total_genes - used_genes
        ug_count.append(not_used)
    return (dict(zip(ug_id, ug_count))) # return dict 

# iteration through a dict's items in random orders, return order  
def iterate_randomly(diction):  
    strainOrder = []        # initialize
    items=list(diction.items()) # create a list of dict items for shuffle
    random.shuffle(items)       # shuffle the tuples
    for strain, identifier in items:
        strainOrder.append(strain)  # save strain 
    return(strainOrder)

# collect the output data 
def gene_contribution(table, strainOrder, ug_dict, base_dict):
    ug = 0              # initialize
    ug_results =[]
    gc_results = []
    table['GC'] = False   # initialize new table column with only False
    for i in strainOrder: # strain order is given shuffled
        # for each strain save True in every row with gene, count Trues
        table.loc[table[i] != '*', 'GC'] = True             
        gc = (len(table['GC'][table['GC'] == True].index)) 
        gc_results.append(gc)    # save gene contribution for output
        ug += ug_dict.get((base_dict.get(i))) 
        ug_results.append(ug) # add unique genes of strain and save
    # return clustered gene contribution and unique genes
    return (gc_results, ug_results, strainOrder)

# set a seed and repeat random iterations of result calculation for a given number of times
def a_few_iterations(iterations, table, table_desc, base, ug_dict):
    random.seed(1234)   # a basic seed  
    data_gc = []        # initialize
    data_ug = []    
    save_order = []
    for i in range(iterations):   # for given number of repetitions
        new_order = iterate_randomly(base) # save the current order
        # collect data of current order
        results = gene_contribution(table, new_order, ug_dict, base)       
        data_gc.append(results[0])   # save results
        data_ug.append(results[1])
        save_order.append(results[2])
    return(data_gc, data_ug, save_order)

### 2. step: visualize box-plots
def visualize_box_plots(listOfRep, table, table_description, base, unique_dict, plotfile_name):

    # create figure with axes for 2 plots, set figure-title
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize= (10,5))
    fig.suptitle("Genome-wise gene contribution", fontsize=18, y=0.97)

    # for iteration number in with number of repitions
    for i in range(len(listOfRep)):

        # calculate data
        data_gc, data_ug, save_order = a_few_iterations(listOfRep[i], table, table_description, base, unique_dict)

        # sum gene contribution and unique genes item wise for total
        total = []                      
        for y in range(listOfRep[i]):
            total.append([sum(x) for x in zip(data_gc[y], data_ug[y])])

        # create DataFrame for plotting
        number_genomes = (table.shape[1] - 3) #GC is additional
        df2 = pd.DataFrame(total, columns=list(range(1,number_genomes)))

        # create boxplot
        axes[i] = df2.boxplot(ax=axes[i])

        # create a style for the display of repetition numbers
        box_style = dict(facecolor='white', edgecolor='grey', boxstyle='round')

        # for all axes use the following layout settings
        #for ax in axes:

        axes[i].set_ylim(2900, 6700)             # frame the y-axis
        axes[i].set_xlabel('Number of genomes')  # label the x-axis

        # create the box of repetition numbers and use style sheet
        axes[i].text(0.03, 0.97, 'Number of repetitions: ' + str(listOfRep[i]), transform=axes[i].transAxes, fontsize=10,
                    verticalalignment='top', bbox=box_style)

    # only label first column plot y-axes  
    axes[0].set_ylabel('Gene cluster count')

    # more layout configuration: distances between plots
    plt.subplots_adjust(hspace=0.35, wspace=0.2)
    fig.tight_layout(pad=3.0,rect=[0, 0.03, 1, 0.95])
    plt.savefig(plotfile_name)
    #plt.show()
