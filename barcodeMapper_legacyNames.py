"""Component of BABUSHKA pipeline responsible for mapping barcodes to known-libraries
   
   Note, this particular version uses 'legacy' names (such as Anc80BC_1317) instead
   of current naming conventions

Copyright (c) 2020 Eric Zinn

This code is free software;  you can redstribute it and/or modify it under the terms of the AGPL
license (see LICENSE file included in the distribution)
"""


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
        data = [(int(reads), str(sequence)) for reads, sequence in csv.reader(inFile)]

    return data


def readToBinaryBarcode(sequence, configFile):
    """Take a pruned, trimmed read and return a binary barcode in the format
    [B1, B2,...,BN] according to the dictionary passed within the configFile"""

    codedBarcode = []

    for barcode in configFile['barcodeArray']:
        index = configFile['barcodeIndices'][barcode]
        codedBarcode.append(configFile['barcodeDict'].get(barcode + "_" + sequence[index:index + 3]))

    return codedBarcode


def nucleotidesToInteger(readCounts, configFile):

    codedCounts = []
    integerCounts = []
    listOfIntegers = []

    for entry in readCounts:
        codedCounts.append((entry[0], readToBinaryBarcode(entry[1], configFile)))

    for entry in codedCounts:
        # If the entry doesn't have a valid Binary Barcode (e.g. mismatch)
        # then don't add it to the list of counts
        if None in entry[1]:
            continue
        # Otherwise, add it
        else:
            integerCounts.append((entry[0], configFile['libName'] + str(binaryToB10(entry[1]))))
            listOfIntegers.append(binaryToB10(entry[1]))

    # Time to clean up counts.  The following block of code combines
    # duplicate counts.  Duplicate counts of course occur when there
    # are polymorphisms in the sequence in between barcode sites...
    # ...in the future, we may want to include an option to discard
    # imperfect barcodes...

    combinedCounts = []

    # for x in range(0,configFile['libSize']):
    for x in set(listOfIntegers):
        duplicateList = [entry for entry in integerCounts if entry[1] == configFile['libName'] + str(x)]
        sumOfDuplicates = sum([pair[0] for pair in duplicateList])
        combinedCounts.append((configFile['libName'] + str(x), sumOfDuplicates))

    return combinedCounts


def binaryToB10(binaryList):
    if(-1 in binaryList):
        # Need to handle negative numbers for the case of Anc83...
        nonNegativeList = [0 if x == -1 else x for x in binaryList]
        #print(-int("".join(str(x) for x in nonNegativeList), 2))
        return -int("".join(str(x) for x in nonNegativeList), 2)
    else:        
        return int("".join(str(x) for x in binaryList), 2)


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
parser.add_argument("-x", "--crane_hack", action = "store_true",
                    help = "Include this parameter if you want to use the 'crane hack' to render the ouput compatible with crane's visualizaiton tool")
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
    print("No libraries found in " + libDirectory)

if not controlList:
    print("No controls found in " + controlList)


# Main Loops

for yamlFile in libList:
    with open(libDirectory + yamlFile, 'r') as ymlfile:
        cfg = yaml.load(ymlfile)
        if args['verbose']:
            print("Mapping lib from " + yamlFile)
        barcodeCounts = nucleotidesToInteger(readCounts, cfg)

        if args['crane_hack']:
            # To increase compatibility with Crane's old library, we're going to write a quick
            # hack to include barcodes with '0' counts in the file and organize the barcodes
            # by their numeric representation (in ascending order)
            hackedCounts = []
            barcodeDict = dict(barodeCounts)

            for x in range(0, yamlFile['libSize']):
                # Check to see if a given barcode has been mapped (by referencing a dict)
                if barcodeDict.get(yamlFile['libName'] + str(x)):
                    # If so, append a tuple containing the barcode name and the number of counts to the
                    # 'hackedCounts' list
                    hackedCounts.append((yamlFile['libName'] + str(x), barcodeDict[(yamlFile['libName'] + str(x))]))
                else:
                    # If not, then append a count of 0
                    hackedCounts.append((yamlFile['libName'] + str(x), 0))

            writeCountsToCSV(outputPath + inputPrefix + cfg['libName'] + '.csv', hackedCounts)

        else:
            writeCountsToCSV(outputPath + inputPrefix + cfg['libName'] + '.csv', barcodeCounts)


controlBarcodeCounts = []

for yamlFile in controlList:
    with open(controlDirectory + yamlFile, 'r') as ymlfile:
        if args['verbose']:
            print("Mapping controls from: " + yamlFile)
        cfg = yaml.load(ymlfile)
        controlBarcodeCounts = controlBarcodeCounts + countControls(cfg, readCounts)

writeCountsToCSV(outputPath + inputPrefix + 'ControlBC.csv', controlBarcodeCounts)
