import sys 
import os
import argparse
from cli_args import number_str
from subjob import write_subjob
#I want to create a script that ask for the step size and creates a subjob that uses it

#Lets assume we are on the directory where i can store all the data 

# Lets do a rework here
import json
from jinja2 import Environment, FileSystemLoader
import numpy as np
# We need jinja and json here


parser = argparse.ArgumentParser(description="Loads the shell tests directories and saves the data in a single file")

parser.add_argument("--folderpath", type=str, required=True, help="path to the folder where the shell tests are stored")
parser.add_argument("--p", type=int, required=True, help="power of the shell")
parser.add_argument("--shell_width", type= float, required=True, help="width of the shell"
                    )
args = parser.parse_args()

folderpath = args.folderpath
p = args.p
shell_width =args.shell_width

def main():
    main_file = "../Results/Shell_test/CoverageSum.txt"
    write_main = open(main_file,"a+")
    f = open(folderpath +"/Coverage_data.txt","r+")
    f.readline()
    # OK so now we can read the data 
    line = f.readline()
    while line:
        write_main.write("{} {} {} \n".format(p,shell_width,line))
        line = f.readline()
    write_main.close()
    f.close()


main()