"""Component of BABUSHKA pipeline responsible for counting barcodes

Copyright (c) 2020 Eric Zinn

This code is free software;  you can redstribute it and/or modify it under the terms of the AGPL
license (see LICENSE file included in the distribution)
"""

import numpy as np
import pandas as pd
import csv
import argparse
from Bio import SeqIO


def loadFastQ(pathToFastq):
    allReads = SeqIO.parse(pathToFastq, "fastq")
    return allReads


def trim_adaptors(record, upstreamAdaptor, downstreamAdaptor, min_len):
    """ Trims perfector adaptor sequences, checks read length."""

    len_up_adaptor = len(upstreamAdaptor)  # cache this for later
    len_record = len(record)  # cache this for later)

    indexUp = record.seq.find(upstreamAdaptor)
    indexDown = record.seq.find(downstreamAdaptor)

    if indexUp == -1 or indexDown == -1:
        # One or more adaptors not found, so don't keep it
        return ""

    if indexDown - indexUp - len_up_adaptor > min_len:
        # after trimming this will still be long enough
        trimmedSequence = record[indexUp + len_up_adaptor:indexDown]

    else:
        trimmedSequence = ""

    return trimmedSequence


def writeCountsToCSV(pathToFile, readCounts):

    with open(pathToFile, 'w', newline = '') as outFile:
        writer = csv.writer(outFile)
        for entry in readCounts:
            writer.writerow(entry)

# Begin agument parser options


parser = argparse.ArgumentParser()

parser.add_argument("-p", "--phred_threshold", type = int,
                    help = "minimum phred quality (default: 30)", default = 30)
parser.add_argument("-r", "--read_threshold", type = int,
                    help = "minimum read length (default: 0)", default = 0)
parser.add_argument("-o", "--output_file", type = str,
                    help = "location of the output file (default: seqCounts.csv)", default = "seqCounts.csv")
parser.add_argument("-m", "--maximum_reads", type = int,
                    help = "Maximum number of reads to process before terminating")
parser.add_argument("-v", "--verbose", action = "store_true",
                    help = "Set verbose (more text. Useful for logging/debugging/progress monitoring")
parser.add_argument("-u", "--upstream_adaptor", type = str,
                    help = "Upstream adaptor sequence to search for", default = "AAGCTT")
parser.add_argument("-d", "--downstream_adaptor", type = str,
                    help = "Downstream adaptor sequence to search for", default = "GCGGCCGC")

requiredNamed = parser.add_argument_group("Required Named Arguments")
requiredNamed.add_argument("-i", "--input_file", type = str,
                           help = "location of the input file (fastq only at the moment)")


args = vars(parser.parse_args())
# print(args)

# Check to make sure the user supplied an input file
if not args['input_file']:
    print('Please provide a valid fastq file')
    quit()


print("Loading File: " + args['input_file'])
sequenceGenerator = loadFastQ(args['input_file'])   

#upstreamAdapter = "AAGCTT"
upstreamAdapter = args['upstream_adaptor']
#downstreamAdapter = "GCGGCCGC"
downstreamAdapter = args['downstream_adaptor']
readLengthThreshold = args['read_threshold']
qualityScoreFilter = args['phred_threshold']

uniqueList = []
countList = []
readCounter = 0
missingAdapterCount = 0
qualityFilterCount = 0

for seqRecord in sequenceGenerator:
    # Print progress to council every 100,000 reads
    if args['verbose']:
        if readCounter % 100000 == 0:
            print("Processed %i reads" % readCounter)

    if args['maximum_reads']:
        if readCounter >= args['maximum_reads']:
            break
    # First, trim the read down to size
    trimmedSeq = trim_adaptors(seqRecord, upstreamAdapter, downstreamAdapter, readLengthThreshold)
    if len(trimmedSeq) == 0:
        readCounter += 1
        missingAdapterCount += 1
        continue

    # Next, test to see if it's of sufficient quality
    if min(trimmedSeq.letter_annotations["phred_quality"]) < qualityScoreFilter:
        readCounter += 1
        qualityFilterCount += 1
        continue

    sequenceString = str(trimmedSeq.seq)

    if sequenceString not in uniqueList:
        uniqueList.append(sequenceString)
        countList.append(1)

    else:
        countList[uniqueList.index(sequenceString)] += 1

    readCounter += 1

print('Done! Processed a total of %i reads' % readCounter)
print('There were %i reads which were missing one or more adapters, or were too small' % missingAdapterCount)
print('There were %i reads which were filtered out because they had a Q score less than the threshold' % qualityFilterCount)
seqCounts = list(zip(countList, uniqueList))
outputFilepath = args['output_file']
print("Writing file to: " + outputFilepath)
writeCountsToCSV(outputFilepath, seqCounts)
