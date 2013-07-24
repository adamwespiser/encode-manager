#!/share/bin/python 
import os, sys, fileinput, pprint, subprocess, optparse
from optparse import OptionParser


home = "/home/wespisea/"

lncDir="/home/wespisea/work/research/researchProjects/encode/data/"

bigBedToBed = "/home/wespisea/bin/bedtools-2.17.0/bin/bigBedToBed"
intersectBed = "/home/wespisea/bin/bedtools-2.17.0/bin/intersectBed"
outdir = "/home/wespisea/scratch/encodePeakSeq-lncRegions/"

parser = OptionParser()
parser.add_option("-i", "--inputFile", dest="inFile", help="read FILE", metavar="INFILE")
parser.add_option("-o", "--outFile", dest="outFile", help="read FILE", metavar="INFILE")
parser.add_option("-l", "--lncList", dest="lncFile", help="list of lncs FILE", metavar="LNCFILE")
(options, args) = parser.parse_args()

def clusterRun(command, mem=4, cores=2):
	cluster = "runJob -m %s -c %s -i \"%s\"" % (mem, cores, command )
	#os.system(cluster)

def swapFileEnding(file,newEnding):
	(base,tail) =  os.path.split(file)
	oldEnding = ".%s" % (tail.split(".")[-1])
	newFile = "/".join([base, tail.replace(oldEnding,newEnding)])
	return newFile

def createBedSummary(bedFile,lncFile,debug=False):
	(base,tail) = os.path.split(bedFile)
	target1 = "/".join([base, tail.replace(".bed",".lnclist"  )])
	target2 = "/".join([base, tail.replace(".bed",".summary"  )])
	# "cat ~/all_lncHits.tab ~/tmp_lncHits.tab |sort| uniq -c | sed 's/     //g' | awk  '{a=$2-1;print $2, "\t", $1-1}' | sort > tmp"	
	#target1 = "~/scratch/tmp"
	cmd1 = "cut -f 14 %s > %s" % (bedFile, target1)
	cmd2 = "cat %s %s | sort | uniq -c | sed \'s/  //g\' | awk \'{print \$2, \\\"\\\\t\\\", \$1-1}\'| sort | sed \'s/ //g\' > %s " % (lncFile,target1,target2)
	cmdAll= "%s && %s" %  (cmd1, cmd2)
	test1 = "wc -l %s |cut -f 1 -d \' \' && cat %s |sed \'s/  //g\'|cut -f 2 | sort | uniq -c| perl -e \'$sum=0;while(<>){@a = split / +/, $_;$sum += ($a[1] * $a[2]),\"\\n\";}print $sum.\"\\n\";\'|cut -d \" \" -f 1" % (target1,target2)
	#teststr= "test 1:\na=`wc -l %s |cut -f 1 -d \" \"` && b=`cat %s %s | sort | uniq -c | sed \'s/  //g\' | awk \'{print $2, \"\\t\", $1-1}\'| cut -d \" \" -f 3|sort|uniq -c| perl -e \'$sum=0;while(<>){@a = split / +/, $_;$sum += ($a[1] * $a[2]),\"\n\";}print $sum.\"\\n\";\'|cut -d \" \" -f 1`" % (target1,lncFile, target1)
	print "echo ",tail," && ", test1
	clusterRun(cmdAll)
 	print "\n" 

if __name__ == "__main__":
	inFile  = options.inFile
	outFile = options.outFile
	lncListFile = options.lncFile

	fileOut = open(outFile, 'w')
	inread = fileinput.input(inFile, mode="r")
	fileArray = [line for line in inread]
	header = fileArray[0].strip()
	fileArray.pop(0)
	#fileOut.write(header)

	lncs = ["body","promoter_distal","promoter_proximal","all_regions"]
	foundSum = 0
	totalSum = 0
	found = 0
	total = 0
	for line in fileArray:
		total = total + 1
		line = line.strip("\n")
		(filename,dccAccession,size,cell,antibody,dataType,view,filename,localfile,geoSampleAccession,md5sum,dateSubmitted,lab,dataVersion,url,softwareVersion,exprdir,controlId,type,dateResubmitted,grant,setType,control,labExpId,treatment,protocol,lncRegion,crossFile,crossFileLoc) = line.split('\t')	

		lncFile = "%slncRNA_%s.tab" % (lncDir,lncRegion)
		targetFile = "%s_%s_%s_%s-X-lnc_%s.bed" % (dccAccession,cell,antibody,lab,lncRegion)
		target = "%s%s" % (outdir,targetFile)
		target = target.replace("(", "\(")
		target = target.replace(")", "\)")
		command = "%s %s stdout | %s -wo -f 0.5 -a stdin -b %s > %s" % (bigBedToBed,localfile,intersectBed,lncFile,target)
		cluster = "runJob -m 5 -c 2 -i \"%s\"" % command

		totalSum = totalSum + 1
		summary=swapFileEnding(crossFileLoc,".summary")
		print summary
		if os.path.isfile(summary):
			foundSum = foundSum + 1
		next
		if os.path.isfile(crossFileLoc):
			found = found + 1
		else:
			print "run job: ",cluster
			print "region: ",lncRegion
			print crossFileLoc
			if os.path.isfile(localfile):
				os.system(cluster)
		createBedSummary(crossFileLoc,lncListFile)	
		
		#outline = "\t".join([line,lncRegion,targetFile,"".join([target,"\n"])    ])
		#fileOut.write(outline)
		
	print "found ", foundSum, "out of ",totalSum,"summary files from last run\n\n"
	print "found ", found, "out of " , total, " total expirements"
