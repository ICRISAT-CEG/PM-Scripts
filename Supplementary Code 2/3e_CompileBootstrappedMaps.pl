#! /usr/bin/perl -w

#Compile results of bootstrapped maps
#Arguments: (0+) Input files, (1) Output file (=ragged arrays since some markers may get dropped from some analyses)

use strict;
my ($snpID, $cmID) = (0,1);	#Columns with snp ID and cM position

my @ins=@ARGV;
my $out = pop @ins;

my %snps;
print "Processing " . scalar(@ins) . " input files:\n";
foreach my $in (@ins){
	print "\t$in\n";
	open(IN, $in) or die "Cannot open input file $in\n";
	my $line;
	#Find start of linkage group
	do{
		$line = <IN>;
	}until ($line =~ m/;BEGINOFGROUP/ || eof);	#Read until hit beginning of group or end of file
	
	#Go through linkage group and load into temporary hash
	my %temp;
	my $max=0;
	while(my $line = <IN>){
		if($line =~ m/;ENDOFGROUP/){last;}	#Break when hit end of group
		chomp $line;
		my ($snp, $cm) = (split/\t/, $line)[$snpID, $cmID];
		#push @{$temp{$snp}}, $cm;
		$temp{$snp} = $cm;
		if($cm > $max) {$max = $cm;}	#Readjust to have
		#print "Max=$max\n";
	}
	close IN;
	
	#Process and save all results
	#print "Max = $max\n";
	foreach my $snp (keys %temp){
		#push @{$snps{$snp}}, ($temp{$snp} / $max);
		#print "Pushing " . $temp{$snp} . "/$max = " . ($temp{$snp} / $max) . "\n";
		$snps{$snp} -> {$in} = ($temp{$snp} / $max)
	}
	#if(scalar(keys %temp) != scalar(keys %snps)){
	#	print "\t\tWARNING!!! Number of SNPs in this file does not match master list. Will cause downstream errors.\n";
	#	
	#}
}

print "Outputting results to $out\n";
open(OUT, ">$out") or die "Cannot open output file $out\n";
foreach my $snp (keys %snps){
	#print OUT "$snp\t" . join("\t", @{$snps{$snp}}) . "\n";
	print OUT "$snp";
	foreach my $in (@ins){
		if(exists $snps{$snp} -> {$in} ){
			print OUT "\t" . $snps{$snp} -> {$in};
		}else{
			print OUT "\tNA";
		}
	}
	print OUT "\n";
}
close OUT;

exit;