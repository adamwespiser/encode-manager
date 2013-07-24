#!/usr/bin/perl
use strict;
use warnings;
#convert lncRNA list from derrien 2012, Genome Research paper to BED format by finding the lnc Coords in genecode V10.

my $pp = 2000; 
my $pd = 10000;

my $bodyFile = "../data/lncRNA_body.bed"; 
my $ppFile = "../data/lncRNA_promoter_proximal.bed"; 
my $pdFile = "../data/lncRNA_promoter_distal.bed"; 
my $allFile = "../data/lncRNA_all_regions.bed";


open my $bodyOut, '>', $bodyFile; 
open my $ppOut, '>', $ppFile; 
open my $pdOut, '>', $pdFile; 
open my $allOut, '>', $allFile; 
open my $lncConf, '>', "../data/lncRNA_List.txt";

my $gencodeGTF = "/home/wespisea/scratch/gencode/Homo_sapiens.GRCh37.70.gtf";
my $gencodeLncGTF = "/home/wespisea/scratch/gencode/gencode.v7.long_noncoding_RNAs.gtf";
my $lncRNAFile = "/home/wespisea/work/research/researchProjects/encode/data/Gencode_lncRNAsv7_summaryTable_05_02_2012.tab"; 

open my $in, '<', $lncRNAFile;
my $head = <$in>;
while(defined(my $line = <$in> )){
	$line =~ s/\n//g; 
	my @line = split "\t", $line; 
	my $trans = $line[1];
	my @transArray= split '\.' , $trans;
#	print $transArray[1]."<  >".$transArray[0]."\n";
	my $cmd = `grep $transArray[0] $gencodeLncGTF|head -n 1|cut -f 7`;
	$cmd =~ s/\n//g; 
	my $dir = $cmd;	
	my $chr = (split ':', $line[2])[0];
	my ($start,  $stop) = split '-', (split ':', $line[2])[1];


	my $startpp, my $startpd;
	my $stoppp,  my $stoppd;
	my $startall,  my $stopall;
	if ($dir ~~ '+'){
		$startpp = $start - $pp;
		$stoppp  = $start; 	
		$startpd = $start - $pd;
		$stoppd  = $start - $pp; 
		$stopall = $stop;
		$startall = $startpd;
	} 
	elsif ($dir  ~~ '-'){

		$startpp = $stop;
		$stoppp  = $stop + $pp; 	
		$startpd = $stop + $pp;
		$stoppd  = $stop + $pd; 
		$startall = $start;
		$stopall = $stop + $pd;
	}
	else {print "notFound: $line\n";
		my $check  = "grep $transArray[0] $gencodeLncGTF";
		my $result = `grep $transArray[0] $gencodeLncGTF`;
		print "$check = $result\n\n\n";      
	
	
	
		next; }
####### REMOVE 
	print $lncConf $trans."\n";
	next;
###### END REMOVE
	print $bodyOut $chr."\t".$start."\t".$stop."\t".$trans."\t"."0\t".$dir."\n"; 
	print $ppOut $chr."\t".$startpp."\t".$stoppp."\t".$trans."\t"."0\t".$dir."\n"; 
	print $pdOut  $chr."\t".$startpd."\t".$stoppd."\t".$trans."\t"."0\t".$dir."\n"; 
	print $allOut  $chr."\t".$startpd."\t".$stoppd."\t".$trans."\t"."0\t".$dir."\n"; 



}
