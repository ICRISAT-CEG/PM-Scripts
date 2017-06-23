#! /usr/bin/perl -w

#Convert an Rqtl "csvr" (rotated CSV) file back into a Hapmap (with dummy codings)
#Arguments: (0) Rqtl csvr file, (1) Hapmap output file

use strict;
use POSIX;

my ($in, $out) = @ARGV;
my ($snpID, $chrID, $posID, $startID) = (0..3);	#Data columns

print "Converting $in to hapmap format\n";
open(IN, $in) or die "Cannot open input file $in\n";
open(OUT, ">$out") or die "Cannot open output file $out\n";

#Process header
my $header=<IN>;
my @headdata=split /[,\t]/, $header;	#Split on either comma or tab, depending on how formatted
print OUT "rs#	alleles	chrom	pos	strand	assembly#	center	protLSID	assayLSID	panelLSID	QCcode	" . join("\t", @headdata[$startID..(scalar @headdata -1)]);

#Process rest of text
my $filler="	+	NA	NA	NA	NA	NA	NA	";
while(my $line = <IN>){
	chomp $line;
	#Change allele codings
	$line =~ s/AA/A/g;	#Homozygous A
	$line =~ s/BB/C/g;	#homozygous B
	$line =~ s/B/C/g;	#homozygous B (single-letter)
	$line =~ s/AC/M/g;	#Het 1
	$line =~ s/AB/M/g;	#Het 1
	$line =~ s/CA/M/g;	#Het 2
	$line =~ s/BA/M/g;	#Het 2
	$line =~ s/H/M/g;	#Het (single letter)
	$line =~ s/-/N/g;	#Missing
	my @data = split /[,\t]/, $line;
	
	print OUT "$data[$snpID]\tA/C\t$data[$chrID]\t" . floor($data[$posID] * 1e5) . $filler . join("\t", @data[$startID..(scalar @data -1)]) . "\n";
}

close IN;
close OUT;


exit;