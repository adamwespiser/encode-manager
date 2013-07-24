#!/share/bin/python

import os, sys, fileinput, pprint, subprocess,optparse
from optparse import OptionParser

home = "/home/wespisea/"
lncDir="/home/wespisea/work/research/researchProjects/encode/data/"
parser = OptionParser()
parser.add_option("-i", "--inputFile", dest="inFile", help="read FILE", metavar="INFILE")
parser.add_option("-o", "--outputFile", dest="outFile", help="output to this FILE", metavar="OUTFILE")
(options, args) = parser.parse_args()

#usual input -> /home/wespisea/work/research/researchProjects/encode/data/gencodeV7lncRNAexpr_stats.tab
# "filename NHEK K562 HELAS3 GM12878 H1HESC"
global cellTypes 
cellTypes = ["NHEK", "K562", "HELAS3", "GM12878", "H1HESC"]

global cellTypeHash
cellTypeHash = {0 : "NHEK", 1 : "K562", 2: "HELAS3", 3 : "GM12868", 4: "H1HESC" }

def processLine(line):
	arr = line.split('\t')
	arr[5] = arr[5].strip()
	for i in xrange(1,6):
		try:
			arr[i] = float(arr[i])
		except ValueError:
			arr[i] = arr[i]
	return arr

def printMat(mat):
	#print mat	
	for row in mat:
		print "\t".join(map(lambda x: str(x),row))


def applyFnsToArray(array):
	index = len(array[0])
	array[index] = label
	for line in array[1:]:
		argArray = array[1:6]
		result = fn

if __name__ == "__main__":
	inFile  = options.inFile
	outFile = options.outFile
	#fileOut = open(outFile, 'w')
	inread = fileinput.input(inFile, mode="r")
	fileArray = [processLine(line) for line in inread]
	header = fileArray[0]
	#fileOut.write(header)
	#print "\t".join(header)
	printMat(fileArray)

