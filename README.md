# proteinortho_pangenome
A command line tool for a pangenome analysis from proteinotho calculated cluster table and description files using python 3.8. Add gene descriptions to the tabular output of cluster tables from proteinortho calculations. Retrieve a tabular pangenome analysis of various proteinortho cluster tables of from one calculation. So now you can compare the various cluster output tables from proteinorto.  



# Usage 

usage: main.py [-h] [-i INPUTFILE [INPUTFILE ...]] [-d DESCRIPTION] [-N] [-A] [-V] [-G] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFILE [INPUTFILE ...], --inputfile INPUTFILE [INPUTFILE ...]
                        path to tabular proteinortho output file(s) required for options -A/-G/-N/-V
  -d DESCRIPTION, --description DESCRIPTION
                        path to tabular proteinortho description file
  -N, --NAMES           add gene names to IDs for the tabular input file
  -A, --ANALYSE         tabular pangenome analysis
  -V, --VISUALIZE       visualization of the pangenome composition as pie or bar plot
  -G, --GENE_CONTRIBUTION
                        box plot of gene contribution with each genome
  -v, --version         show program version


# Examples

