#! /use/bin/perl -w

#Convert an MSTmap input and output file pair into an Rqtl file
#Arguments: (0) Input hapmap (in ABH format), (2) Output file name for R/qtl

use strict;

my ($in, $out) = @ARGV;
my ($snpID, $chromID, $posID, $startcol)=(0,2,3,11);	#Columns whith needed information

#Set up arrays to hold everything
#Load input marker set
print "Processing $in to R/qtl csvr format (tab-delimited)\n";
open(IN, $in) or die "Cannot open input file $in\n";
open(OUT, ">$out") or die "Cannot open output file $out\n";
##Process header
my $header=<IN>;
chomp $header;
my @headdata=split/\t/, $header;
print OUT "taxon\t\t\t" . join("\t", @headdata[$startcol..(scalar @headdata - 1)]) . "\n";
while(my $line=<IN>){
	chomp $line;
	my @data = split /\t/, $line;
	print OUT "$data[$snpID]\t$data[$chromID]\t$data[$posID]\t" . join("\t", @data[$startcol..(scalar @data - 1)])  . "\n";
}
close IN;
close OUT;

exit;