#! /usr/bin/perl -w

#Apply ad-hoc filters to VCF file (sort of a cheap workaround, but VCFtools doesn't have a maxGC filter option)
#Note: Filter options currently hardocded
#Arguments: (0) Input file, (1) Output file, (2) Comma-separated list of acceptable GC scores

use strict;

my ($in, $out, $quals) = @ARGV;

my $startcol=9;	#Column where genotypes start
my ($depthID,$qualID) = (2,3);	#Positions in VCF annotation with quality and depth information
#Filter options
my ($minQual, $maxQual) = (89, 99);
my $missing="\\.\\/\\.";	#Missing code, for regex

#Process acceptable quality scores
my %goodqual;
my @qualities= split/,/, $quals;
foreach my $q (@qualities){
	$goodqual{$q}=1;
}

print "Processign $in\n";
open(IN, $in) or die "Cannot open input file $in\n";
open(OUT, ">$out") or die "Cannot open output file $out\n";
my ($counter,$changed,$total)=(0,0,0);
while(my $line=<IN>){
	$counter++;
	if($counter % 10000==0){
		print "\tProcessed $counter lines\n";
	}
	#if($counter>1000){last;}	#For debugging
	if($line =~ m/^#/){	#Header rows
		print OUT $line;
		next;
	}
	chomp $line;
	my @data = split /\t/, $line;
	for(my $i=$startcol; $i<scalar @data; $i++){
		$total+=1;
		#print "Checking $data[$i] for $missing\n";
		if($data[$i] =~ m/^$missing/){next;	}	#Skip missing data (either missing entirely, or set to missing by VCFtools)
		#if($data[$i] =~ m/^\.\/\./){next;	}	#Skip missing data (either missing entirely, or set to missing by VCFtools)
		#print "\tProcessing $data[$i]\n";
		my ($depth, $qual) = (split /:/, $data[$i])[$depthID, $qualID];
		#print "$data[$i] has depth $depth and quality $qual\n";	#For debugging
		if(!exists $goodqual{$qual}){	#Set unacceptable GCs to missing
			$data[$i] = $missing;
			$changed+=1;
		}
	}
	print OUT join("\t", @data) . "\n";
}
print "$changed sites removed out of $total total (portion=" . ($changed/$total) . ")\n";

close IN;
close OUT;

exit;