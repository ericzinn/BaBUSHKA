import numpy as np
import pandas as pd
import csv
import argparse
import yaml

from Bio import SeqIO
from collections import namedtuple
from collections import Counter
from os import listdir
from os.path import isfile, join, basename

def readCSVtoCounts(pathToFile):
    with open(pathToFile) as inFile:
        data=[(int(reads),str(sequence)) for reads, sequence in csv.reader(inFile)]

    return data

def readToBinaryBarcode(sequence, configFile):
    """Take a pruned, trimmed read and return a binary barcode in the format
    [B1, B2,...,BN] according to the dictionary passed within the configFile"""

    codedBarcode = []

    for barcode in configFile['barcodeArray']:
        index = configFile['barcodeIndices'][barcode]
        codedBarcode.append(configFile['barcodeDict'].get(barcode+"_"+sequence[index:index+3]))
    
    return codedBarcode

def nucleotidesToInteger(readCounts, configFile):
    
    codedCounts = []
    integerCounts = []

    for entry in readCounts:
        codedCounts.append((entry[0],readToBinaryBarcode(entry[1], configFile)))
    
    for entry in codedCounts:
        if None in entry[1]:
            continue
        else:
            integerCounts.append((entry[0],configFile['libName']+str(binaryToB10(entry[1]))))
            
    # Time to clean up counts.  The following block of code combines
    # duplicate counts.  Duplicate counts of course occur when there
    # are polymorphisms in the sequence in between barcode sites...
    # ...in the future, we may want to include an option to discard
    # imperfect barcodes...
    
    combinedCounts = []
    
    for x in range(0,configFile['libSize']):
        duplicateList = [entry for entry in integerCounts if entry[1] == configFile['libName']+str(x)]
        sumOfDuplicates = sum([pair[0] for pair in duplicateList])
        combinedCounts.append((configFile['libName']+str(x), sumOfDuplicates))
    
        
    return combinedCounts

def binaryToB10(binaryList):
    return int("".join(str(x) for x in binaryList),2)


def writeCountsToCSV(pathToFile, readCounts):    
    with open(pathToFile, 'w', newline = '') as outFile:
        writer = csv.writer(outFile)
        for entry in readCounts:
            writer.writerow(entry)

def countControls(configFile, readCounts):
    
    countsList = []
    
    for control in configFile['controlBC']:
        try:
            BCCount = [v[0] for v in readCounts if v[1] == configFile['controlBC'][control]][0]
        except:
            BCCount = 0
        countsList.append((control, BCCount))
        
    return countsList

    
    
parser = argparse.ArgumentParser()


parser.add_argument("-v", "--verbose", action = "store_true",
                    help = "Set verbose (more text. Useful for logging/debugging/progress monitoring")
parser.add_argument("-o", "--output_directory", default = "",
                    help = "Path to where you want the output to live")
requiredNamed = parser.add_argument_group("Required Named Arguments")
requiredNamed.add_argument("-i", "--input_file", type = str,
                    help = "location of the input file (csv file of ReadCounts)")
requiredNamed.add_argument("-l", "--lib_directory", type = str,
                    help = "directory containing the yaml Library Files")
requiredNamed.add_argument("-c", "--control_directory", type = str,
                    help = "directory containing the yaml control files")
                    
args = vars(parser.parse_args())

# Begin argument checking.  Future sanitization here

libDirectory = args['lib_directory']
controlDirectory = args['control_directory']
inputFile = args['input_file']
outputPath = args['output_directory']

inputPrefix = basename(inputFile).split('.csv')[0]
if not libDirectory:

    print("Must provide a lib directory")
    quit()

if not controlDirectory:
    print("Must provide a control directory")
    quit()

if not inputFile:
    print("Must provide an input file")
    quit()

    
libList = [f for f in listdir(libDirectory) if f.endswith(".yaml")]
controlList = [f for f in listdir(controlDirectory) if f.endswith(".yaml")]
readCounts = readCSVtoCounts(inputFile)

if not libList:
    print("No libraries found in "+libDirectory)

if not controlList:
    print("No controls found in "+controlList)
    

# Main Loops

for yamlFile in libList:
    with open(libDirectory+yamlFile, 'r') as ymlfile:
        cfg = yaml.load(ymlfile)
        if args['verbose']:
            print("Mapping lib from "+yamlFile)
        barcodeCounts = nucleotidesToInteger(readCounts, cfg)
        writeCountsToCSV(outputPath+inputPrefix+cfg['libName']+'.csv', barcodeCounts)

        
controlBarcodeCounts = []

for yamlFile in controlList:
    with open(controlDirectory+yamlFile, 'r') as ymlfile:
         if args['verbose']:
             print("Mapping controls from: "+yamlFile)
         cfg = yaml.load(ymlfile)
         controlBarcodeCounts = controlBarcodeCounts + countControls(cfg, readCounts)

writeCountsToCSV(outputPath+inputPrefix+'ControlBC.csv', controlBarcodeCounts)