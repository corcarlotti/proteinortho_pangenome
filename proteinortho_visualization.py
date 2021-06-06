# pie chart and bar plot visualization of pangenome 

import pandas as pd
import matplotlib.pyplot as plt
import sys, getopt
import numpy as np

# function for labeling absoulte and relative values
def label_func(pct, array):
    absolute = int(pct/100.*np.sum(array))
    return "{:.0f}\n({:2.2f}%)".format(absolute,pct)

# genefamilies with genes shared by all species (core genome) 
def core_clusters(table):
    #over table collect in array False if row contains '*'
    coreclusters_bool = table.apply(lambda x: False if '*' in list(x)
                                    else True, axis=1) 
    coreclusters_rows = len(coreclusters_bool[coreclusters_bool == True].index) #number of False in array
    return(coreclusters_rows)    

# genefamilies shared by at most four species (accessory genome) 
def acc_clusters(table):
    #over table collect in array True if row contains '*'
    acclusters_bool = table.apply(lambda x: True if '*' in list(x)
                                  else False, axis=1) 
    acclusters_rows = len(acclusters_bool[acclusters_bool == True].index) #number of True in array
    return (acclusters_rows)

# genes that are not included in clusters
def unique_gene_calc(table, table_desc):
    total_genes = len(table_desc['IDs'].unique())
    included_genes = table['Genes'].sum()
    return(total_genes - included_genes )

# build type array for number of species 
def build_type_array(table):
    number_genomes = (table.shape[1] - 3)
    type_array = ["unique genes"]
    for x in range(number_genomes-2):
        type_array += ["accessory gene cluster"]
    type_array += ["core gene cluster"]
    return (type_array)

# match colors to cluster types
def build_color_dict(number_genomes, unique_bool):
    str_for_strain = map(str, (range(1,(number_genomes)+1)))
    colors = ['#1aa3ff','#99d6ff','#007acc']
    if unique_bool: 
        color_array = [colors[0]]
        for x in range(number_genomes-2):
            color_array += [colors[1]]
        color_array += [colors[2]]
    else:
        color_array = [colors[1]]
        for x in range(number_genomes-2):
            color_array += [colors[1]]
        color_array += [colors[2]]
    color_dict= dict(zip(str_for_strain, color_array))
    return (color_dict)

def pangenome_visualize(table, description, plotfilename):

### 1. part: pie charts ###
    
    ## Initializing and building favorized order

    # gather data
    accessory_genes = acc_clusters(table)
    core_genes = core_clusters(table)
    unique_genes = unique_gene_calc(table, description)
    
    # first subplot (3 groups, including unique genes)
    pangenome_parts=['acc', 'core', 'unique']
    pangenome_values=[accessory_genes, core_genes, unique_genes]

    # second subplot (2 groups, unique genes included in accessory)
    pangenome_small_parts=['core', 'acc']
    new_acc= (accessory_genes + unique_genes)
    pangenome_small_values=[core_genes, new_acc]

    ## Create Series, Plots and give Characteristics

    # create main plot figure and axes
    f, axs = plt.subplots(2,2, figsize=(9,7))

    # first subplot
    piechart_series = pd.Series(pangenome_values,
                        index=pangenome_parts)
    axs[0, 0].pie(piechart_series,
        # use autopct of pie() to send values to label_func via labmda function
        autopct=lambda pct: label_func(pct, pangenome_values),
        labels= pangenome_parts,
        colors= ['#1aa3ff','#99d6ff','#007acc'],
        startangle=0)
    axs[0, 0].set_title("core - accessory - unique", fontsize=10)
    axs[0, 0].set_ylabel("")
    axs[0, 0].set_xlabel("cutoff=0.95", fontsize= 7,
                         ha="right", va="baseline")

    # second subplot
    piechart_series2 = pd.Series(pangenome_small_values,
                        index=pangenome_small_parts)
    axs[0, 1].pie(piechart_series2,
        # use autopct of pie() to send values to label_func via labmda function
        autopct=lambda pct: label_func(pct, pangenome_small_values),
        labels= pangenome_small_parts,
        colors= ['#99d6ff','#1aa3ff'],
        startangle=110)
    axs[0, 1].set_title("core - accessory", fontsize=10)
    axs[0, 1].set_ylabel("")
    axs[0, 1].set_xlabel("cutoff=0.95", fontsize= 7,
                         ha="right", va="baseline")

### 2.part: bar plot ###

    # group cluster by number of involved species  
    strain_count = pd.Series((table.groupby(['# Species'],as_index=False).count())['Genes'])

    # create pandas Series object for dataFrame input out of retrieved data
    the_numbers= (pd.Series(unique_genes)).append(strain_count,  ignore_index=True)
   
    ### 2. step visualize
    
    # create dataFrame
    number_genomes= (table.shape[1] - 3)
    barplot_df= pd.DataFrame({'Number of strains':map(str, (range(1,(number_genomes)+1))),
                              'Type':build_type_array(table),
                              'Gene cluster': the_numbers})

    # initialize colors, bool for manner of display
    colors1 = build_color_dict(number_genomes,True) 
    colors2 = build_color_dict(number_genomes,False)
    uni_colors = np.unique(list(colors1.values()))

    # give 1.  and 2. plot data and some characteristics
    barplot_df.plot.bar(ax= axs[1, 0], x='Number of strains', y='Gene cluster', rot= 0, legend= False, color=list(colors1.values()))
    barplot_df.plot.bar(ax= axs[1, 1], x='Number of strains', y='Gene cluster', rot= 0, legend= False, color=list(colors2.values()))

    # for-loop for labels on bars, which tell size of bars
    for i, barplot in enumerate(axs[1]):

        axs[1, i].set_ylabel("Gene cluster count")

        for i, rect in enumerate(barplot.patches):
            height=rect.get_height()
            width=rect.get_width()
            x=rect.get_x()
            y=rect.get_y()
            label_text = int(height)

            label_x= x + width / 2
            label_y= y + height / 2
            barplot.text(label_x, label_y, label_text,
                         ha='center', va='center')

        # hide y-axis numbers and label
        ax1 = plt.axes
        y_axis = barplot.axes.get_yaxis()
        y_axis.set_visible(False)

        barplot.spines["top"].set_visible(False)
        barplot.spines["right"].set_visible(False)
        barplot.spines["left"].set_visible(False)

        barplot.set_xlabel("Strains included in gene cluster")

    # create legend and give legend text and color 1. barplot
    labels1 = barplot_df.Type.unique()
    handles1 = [plt.Rectangle((0,0),1,1, color=uni_colors[i])
                for i , label in enumerate(labels1)]
    axs[1, 0].legend(handles1, labels1)

    # create legend and give legend text and color 2. barplot
    labels2 = barplot_df.Type.unique()[1:]
    handles2 = [plt.Rectangle((0,0),1,1, color=uni_colors[1:][i])
                for i , label in enumerate(labels2)]
    axs[1, 1].legend(handles2, labels2)
    
    # show plot
    f.tight_layout(pad=3.0)
    plt.savefig(plotfilename)
    #plt.show()
