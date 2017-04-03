#! /usr/bin/perl -w

#Purpose: Select a random subset of snps from a set of Hapmap files to generate an MDS plot
#Arguments: (0) Input file; (1) n (total line number to keep); (2) header? (passed as "header"); optional; (3) output file

use strict;

my $infile = shift @ARGV;
my $n = shift @ARGV;
my $outfile = pop @ARGV;
my $doHeader="no";
if(scalar (@ARGV) > 0 && $ARGV[0] eq "header"){
	$doHeader = "yes";
}

#Determine # lines per input file
print "Counting lines in file $infile\n";
my $lines=0;
open(INFILE, $infile) or die "Cannot open input file $infile\n";
while(my $line=<INFILE>) {$lines+=1;}
close INFILE;
if($doHeader eq "yes") {$lines--;}	#don't count the header line
print "\tTotal $lines lines\n";

#Determine lines
print "Determining $n random lines to print\n";
my @randoms;
my %linekey;
if($lines > $n){
	for(my $i=0; $i<$n; $i++){
		$randoms[$i] = int(rand($lines));
		if(exists($linekey{$randoms[$i]})){	#If this line was already chosen, decrement $i and to effectively redo this iteration
			$i--;
		}
		else{
			$linekey{$randoms[$i]}=1;	#Store line position
			#print "Stored line position $randoms[$i]\n";
		}
	}
}else{	# If fewer lines in file than requested
	for(my $i=0; $i<$n; $i++){
		$linekey{$i}=1;
	}
}
#print "\tTotal " . scalar @randoms . " lines to keep out of $n; hash has " . scalar(keys %linekey). " elements\n";

#Output lines to new file
print "Outputting subset\n";
open(INFILE, $infile);
open(OUTFILE, ">$outfile") or die "Cannot open output file $outfile\n";
if($doHeader eq "yes"){
	my $header=<INFILE>;
	print OUTFILE $header;
}
my $linenum=0;
my $total=0;
while(my $line=<INFILE>){
	
	if(exists $linekey{$linenum}){
		print OUTFILE $line;
		$total++;
		$linekey{$linenum}=0;
	}
	$linenum+=1;
}
close INFILE;
close OUTFILE;
print "Total $total lines outputted\n";

if($total < $n){
	foreach my $line (keys %linekey){
		if($linekey{$line}==1){
			print "Warning! Line $line not printed!!\n";
		}
	}
}

exit;