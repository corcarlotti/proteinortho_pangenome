### translate description ID into ID + Description

import sys, getopt

def add_desc(inputfile, description, outputfile):

        # try reading input
        try:
                input_mainfile = open(inputfile, "rt")
                mainfile = input_mainfile.read()     

        #throw an exception if not working
        except IOError:
                print ("Can't read this input file.")
                sys.exit(2) 

        # try reading description
        try:
                input_namelist = open(description, "rt")
                namelist = input_namelist.readlines()
        
        # throw error if not working
        except IOError:
                print ("Can't read this description file.")
                sys.exit(2) 

        # replace ID by ID+Description
        for line in namelist:
                identifier = line.split()[0]
                newline = '%'.join(line.split())+" "
                mainfile = mainfile.replace(identifier,newline)

        # try writing the outputfile
        try:
                input_mainfile = open(outputfile, "w+")
                input_mainfile.write(mainfile)
                input_mainfile.close()
        
        # throw error if not working
        except IOError:
                print ("Please add a name for the output file.")
                sys.exit(2) 
