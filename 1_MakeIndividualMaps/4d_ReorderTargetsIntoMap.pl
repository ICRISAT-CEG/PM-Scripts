#! /usr/bin/perl -w

#Take the Target SNPs and reorder them back into the reference hapmap
#Arguments: (0) Original hapmap input file, (1) Good SNP position file, (2) Bad snps to leave out (3) Output file (hapmap)

use strict;

my ($in, $good, $bad, $out) = @ARGV;
my ($snpID, $chromID, $posID)=(0,2,3);	#Columns in hapmap with chromosome and position number

#Load target SNPs and their best locations
print "Loading best locations of Target SNPs\n";
my %good_snps;
my %neighbor_snps;
open(GOOD, $good) or die "Cannot open best position file $good\n";
<GOOD>;	#Clear header
while(my $line= <GOOD>){
	chomp $line;
	my ($snp, $neighbor, $rsq) = split /\t/, $line;
	if(exists $neighbor_snps{$neighbor}){
		$neighbor_snps{$neighbor} .= "~~~~" . $snp;	#In case more than one SNP at a location, store
	}
	else{
		$neighbor_snps{$neighbor}=$snp;
	}
	$good_snps{$snp}=1;
	#print "Loaded $snp at $site\n";
}
close GOOD;
print "\tLoaded " . scalar(keys %good_snps) . " snps at " . scalar(keys %neighbor_snps)  . " sites\n";

#Go through hapmap file and extract target SNPs to move
my $found=0;
print "Finding target SNPs\n";
open(IN, $in) or die "Cannot open input file $in\n";
while(my $line=<IN>){
	$line =~ m/^(\S+)/;
	my $snp=$1;
	#print "Snp = $snp\n";
	if(exists $good_snps{$snp}){
		$good_snps{$snp} = $line;	#Store information in
		$found++;
	}
}
close IN;
print "\tLoaded $found SNPs total in hapmap to reorder\n";
if($found != scalar(keys %good_snps) ){
	print "\tWARNING! SNPs found do not match best SNPs!!\n";
}
close IN;

#Identify bad SNPs to leave out
open(BAD, $bad) or die "Cannot open file of bad snps $bad\n";
my @bads=<BAD>;
chomp @bads;
close BAD;
my %bad_snps;
foreach my $badsnp (@bads){
	$bad_snps{$badsnp}=1;
}
print "\tLoaded " . scalar(keys %bad_snps) . " bad SNPs to exclude\n";

#Go through and reorder hapmap
print "Reordering hapmap - NOTE: Deleting any markers on chromosome 0 under assumption that they were not mapped\n";
open(IN, $in) or die "Cannot open input file $in\n";
open(OUT, ">$out") or die "Cannot open output file $out\n";
$found=0;
while(my $line = <IN>){
	my ($snp, $chrom, $pos) = (split /\t/, $line)[$snpID, $chromID, $posID];
	if(exists $good_snps{$snp} || exists $bad_snps{$snp} || $chrom eq "0"){
		next;	#skip lines with target SNPs (since supposed to be rearranging them) or bad SNPs
	}
	print OUT $line;	#Print out site, then go on to processing any needed SNPs
	if(exists $neighbor_snps{$snp}){	#If this is where a SNP is supposed to be inserted
		my @snps = split /~~~~/, $neighbor_snps{$snp};
		foreach my $s (@snps){
			my @data = split /\t/, $good_snps{$s};
			$data[$chromID] = $chrom;
			$data[$posID] = $pos;
			print OUT join("\t", @data);
			$found++;
		}
	}
	
}
print "Placed $found SNPs at new locations\n";

close IN;
close OUT;

exit;