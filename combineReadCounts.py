import numpy as np
import pandas as pd
import csv
import argparse
from Bio import SeqIO
from collections import namedtuple
from collections import Counter
from os import listdir
from os.path import isfile, join

def loadReadCounts(pathToFile):
    Row = namedtuple('Row', ('count', 'read'))

    with open(pathToFile, 'r') as f:
        r = csv.reader(f, delimiter = ',')
        readCounts = [Row(*l) for l in r]
    
    return readCounts

def writeCountsToCSV(pathToFile, readCounts):
    with open(pathToFile, 'w', newline = '') as outFile:
        writer = csv.writer(outFile)
        for entry in readCounts:
            writer.writerow(entry)

# Begin agument parser options

parser = argparse.ArgumentParser()

parser.add_argument("-o", "--output_file", default = "totalCounts.csv",
                    help = "The name of your output file (default: totalCounts.csv)")
parser.add_argument("-v", "--verbose", action = "store_true",
                    help = "Set verbose (more text. Useful for logging/debugging/progress monitoring")
requiredNamed = parser.add_argument_group("Required Named Arguments")
requiredNamed.add_argument("-d", "--directory", type = str,
                    help = "location of the csv files to be combined (directory)")

args = vars(parser.parse_args())
dataDirectory = args['directory']

#Check to make sure the user supplied an input file
if not dataDirectory:
    print('Please provide a valid directory')
    quit()
    
fileList = [f for f in listdir(dataDirectory) if f.endswith(".csv")]

if not fileList:
    print("No csv files found in directory specified. Quitting")
    quit()

if args['verbose']:
    print("Files found in directory..."+str(fileList))

frames = []
for file in fileList:
    if args['verbose']:
        print("Loading: "+file)
    frames.append(pd.read_csv(dataDirectory+file, header = None, names = ['count', 'read']))
    
countsDataFrame = pd.concat(frames)
countsDataFrame['Total'] = countsDataFrame.groupby(['read'])['count'].transform('sum')
tempDataFrame = countsDataFrame[['Total','read']].drop_duplicates()

tempDataFrame.to_csv(args['output_file'], index=False, header=False)
print('DONE!')