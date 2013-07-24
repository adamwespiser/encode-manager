#!/usr/bin/perl
#
my $readIn = "/home/wespisea/work/research/researchProjects/encode/data/RnaSeqExpr.tab"; 
my $outfile = "/home/wespisea/scratch/encodeRnaSeq/combinedExpr.tab";
$test = 0;
open my $fh, '<', $readIn;

@cellTypes;
@files;
my $header = <$fh>;
while(defined(my $line = <$fh>)){
 my @arr = split "\t", $line; 
 my $file1 =$arr[10];
 my $cell = $arr[3]; 
 my $rnaExtract = $arr[30];
 next if ($rnaExtract =~ m/longNonPolyA/); # well just get this using the 
 my $file2 = $file1;
 $file2 =~ s/Pap/Pam/g; 
 print $file1."\n";
 print $file2."\n";
 $file1=~s/\.gz//g;
 $file2=~s/\.gz//g;
 my $file11 = $file1."tmp";
 my $file22 = $file2."tmp";
 my $joinFile = $file1;
 $joinFile =~ s/Pap/PamAndPap/g;
 my $cmd1 = "cat $file1 | awk \'{print \$12, \"\\t\",\$6}\'|sed \'s/[\\\";]//g\' > $file11 "; 	
 my $cmd2 = "cat $file2 | awk \'{print \$12, \"\\t\",\$6}\'|sed \'s/[\\\";]//g\' > $file22 "; 	
 my $cmd3 = "join $file11 $file22  | awk \'{print \$1,(\$2+\$3)/2}\' > $joinFile";
 print($cmd1);
 print($cmd2);
  print($cmd3."\n");
 if ($test == 2){
	system($cmd1);
  system($cmd2);
  system($cmd3);
 }
 push(@cellTypes,$cell);
 push(@files, $joinFile);
}
my $ct = join "\t", @cellTypes;
$ct = "transcript_id\t".$ct."\n";
my $fs = join ' ',@files;

print $ct."\n".$fs."\n";
$tmp1 = "~/scratch/tmp1";
$tmp2 = "~/scratch/tmp2";

$cmd = "join $files[0] $files[1] > $tmp1";
$cmd = $cmd." && join $tmp1 $files[2] > $tmp2";
$cmd = $cmd." && join $tmp2 $files[3] > $tmp1";
$cmd = $cmd." && join $tmp1 $files[4] > $outfile";
#print $cmd."\n";

system($cmd) if $test == 1;
print $fs."\n"


