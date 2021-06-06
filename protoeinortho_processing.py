
### functions for proteinortho processing - they all return a value ###

import sys, getopt
import pandas as pd

# total genes = number of distinct identifiers in description list
def total_genes(table_desc):
    totalgenes = len(table_desc['IDs'].unique()) #creates array without duplicates 
    return (totalgenes)

#  clustered genes = number of genes of proteinortho output
def included_genes(table):
    sum_genes = table['Genes'].sum() #sums all genes counted per cluster
    return(sum_genes)

# gene families = number of rows in proteinortho cluster table
def number_clusters(table):
    geneclusters = table.shape[0] #gives back the number of rows in a table
    return(geneclusters)

# genefamilies with genes shared by all species (core genome) 
def core_clusters(table):
    coreclusters_bool = table.apply(lambda x: False if '*' in list(x) else True, axis=1) #Bool-Array: False if the row contains '*'
    coreclusters_rows = len(coreclusters_bool[coreclusters_bool == True].index) #number of True in array
    return(coreclusters_rows)    

# genefamilies with genes shared by at most four species (accessory genome) 
def acc_clusters(table):
    acclusters_bool = table.apply(lambda x: True if '*' in list(x) else False, axis=1) #Bool-Array: True if row contains '*'
    acclusters_rows = len(acclusters_bool[acclusters_bool == True].index) #number of True in array
    return (acclusters_rows)

# gene families not consisting of hypothetical proteins only:  count if #'hypothetical%proteins'in row == #Genes  
def cluster_func(table):
    i = table.columns[3:][0] #select a genome column
    line = pd.Series.to_string(table[table['# Species'] == (table.shape[1] -3)][i][:1], index=False) # select row with id

    table['sum_hp'] = 0
    for column in table.columns[3:]:
        table['sum_hp'] += table[column].astype(str).str.count("hypothetical%protein")
    if "%" in line:
        hp_bool = table.apply(lambda x: False if x['Genes'] == x['sum_hp'] else True, axis=1) #Bool-Array: False if sum_hp == #Genes 
        hp_rows = len(hp_bool[hp_bool == True].index) #number of Trues in array
    else:
        hp_rows = 0 #in case there are no annotated genes
    return(hp_rows)

# total number of core genes counted in core gene clusters
def core_genes(table):
    #print(table.shape[1] - 3)
    coregenes = table.loc[table['# Species'] == (table.shape[1] - 3), 'Genes'].sum() # sum #genes if #species == 5
    return (coregenes) 

# number of accessory genes counted in accessory gene clusters
def acc_genes(table, dif):
    accgenes = table.loc[table['# Species'] < (table.shape[1] - 3), 'Genes'].sum() # sum #genes if #species < 5
    return (accgenes-dif)

# genes that are not included in clusters
def unique_genes(table, table_desc, dif):
    unigenes = (total_genes(table_desc)-included_genes(table))+dif #if there are multiple appearances, subtract number of multiple genes 
    return(unigenes)  

# check for multiple occurrences of a gene identifier, return number of multiple genes
def check_multiple_occurrences(table):

    # check if replaced table or not -> cut value
    i = table.columns[3:][0]        # select a genome colum
    line = pd.Series.to_string(table[table['# Species'] == (table.shape[1] - 3)][i][:1], index=False) # select row with id
    if "%" in line: #if replace identifier included, add additionally included space to separator
        cut = ' ,'
    else:
        cut = ','

    # initialize
    number_of_dupl_genes = 0
    dupl_genes = []

    # for all strains
    for column in (table.columns[3:]):
        
        #drop all values "*", split row at ", " and create list with both columns 
        new_series = pd.Series(table.loc[table[column] != "*"][column].values).reset_index(drop=True)
        expanded = new_series.str.split(cut, expand=True)
        new_list=pd.Series(expanded[0].append(expanded[1]).reset_index(drop=True)).dropna().reset_index(drop=True)

        #extract duplicated values, count unique duplicates
        dupl_genes = new_list[new_list.duplicated() == True]
        number_of_dupl_genes += len(dupl_genes.unique())
        
    return(number_of_dupl_genes) 

### table creation and filling ###

def create_df():

    # prepare index
    process_results = ['Total genes from all input genomes',
                       'Genes in gene clusters',
                       'Gene clusters',
                       'Unique genes',
                       'Core gene clusters',
                       'Accessory gene clusters',
                       'Core genes in clusters',
                       'Accessory genes in clusters',
                       'Gene cluster with annotated function']

    # create empty dataframe
    new_df = pd.DataFrame(index= process_results)
    return (new_df)

    
def fill_df(tab, table_desc, dataframe, tablename):
    # create a new column in dataframe with given column name and values collected in output-array
    column_name = tablename

    # save number of genes with multiple occurrences of gene identifiers in cluster table 
    multiple_genes =check_multiple_occurrences(tab)
    if multiple_genes != 0:
        print("WARNING: There are multiple occurrences of genes in " +tablename+ "\nPlease check your input data or proteinortho calculations.\n" )

    output=[total_genes(table_desc),
            included_genes(tab),
            number_clusters(tab),
            unique_genes(tab, table_desc, multiple_genes),
            core_clusters(tab),
            acc_clusters(tab), 
            core_genes(tab),
            acc_genes(tab, multiple_genes),
            cluster_func(tab)]

    dataframe[column_name] = output
    return(dataframe)
