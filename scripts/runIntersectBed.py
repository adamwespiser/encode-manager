#!/share/bin/python

import os, sys, fileinput, pprint, subprocess,optparse
from optparse import OptionParser

home = "/home/wespisea/"
lncDir="/home/wespisea/work/research/researchProjects/encode/data/"
bigBedToBed = "/home/wespisea/bin/bedtools-2.17.0/bin/bigBedToBed"
intersectBed = "/home/wespisea/bin/bedtools-2.17.0/bin/intersectBed"
outdir = "/home/wespisea/scratch/encodePeakSeq-lncRegions/"
parser = OptionParser()
parser.add_option("-i", "--inputFile", dest="inFile", help="read FILE", metavar="INFILE")
parser.add_option("-o", "--outputFile", dest="outFile", help="output to this FILE", metavar="OUTFILE")
(options, args) = parser.parse_args()


if __name__ == "__main__":
	inFile  = options.inFile
	outFile = options.outFile
	fileOut = open(outFile, 'w')
	inread = fileinput.input(inFile, mode="r")
	fileArray = [line for line in inread]
	header = fileArray[0].strip()
	fileArray.pop(0)
	header = "\t".join([header,"lncRegion","crossFile","crossFileLoc\n"])
	fileOut.write(header)

	lncs = ["body","promoter_distal","promoter_proximal","all_regions"]
	for line in fileArray:
		(filename,dccAccession,size,cell,antibody,dataType,view,filename,localfile,geoSampleAccession,md5sum,dateSubmitted,lab,dataVersion,url,softwareVersion,exprdir,controlId,type,dateResubmitted,grant,setType,control,labExpId,treatment,protocol) = line.split('\t')	
		line = line.strip("\n")
		for i in range(0,3):
			lncRegion = lncs[i]
			lncFile = "%slncRNA_%s.tab" % (lncDir,lncRegion)
			targetFile = "%s_%s_%s_%s-X-lnc_%s.bed" % (dccAccession,cell,antibody,lab,lncRegion)
			target = "%s%s" % (outdir,targetFile)
			command = "%s %s stdout | %s -wo -f 0.5 -a stdin -b %s > %s" % (bigBedToBed,localfile,intersectBed,lncFile,target)
			cluster = "runJob -m 5 -c 2 -i \"%s\"" % command
			#print cluster
			outline = "\t".join([line,lncRegion,targetFile,"".join([target,"\n"])    ])
			fileOut.write(outline)
			if os.path.isfile(targetFile):
				print ".."
			else:
				os.system(cluster)
