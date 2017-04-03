#! /usr/bin/perl -w

#Get stats on the VCF SNPs in terms of genome quality and heterozygosity, since Katie's GC>=98 seems inappropriate for this pop (almost all are hets, meaning almost all are paralogs)
#Arguments: (0) Input file, (1) Output file

use strict;

my($in, $out) = @ARGV;
my $startcol=9;	#Column in VCF file where genotypes start
my ($alleleID, $depthID, $qualID) = (0, 2, 3);#portions of VCF genotypes with needed information


if($in =~ m/.gz$/){
	open(IN, "gunzip $in |") or die "Cannot open zipped input file $in\n";
}else{
	open(IN, $in) or die "Cannot open input file $in\n";
}

open(OUT, ">$out" ) or die "Cannot open output file $out\n";
print OUT "GC\tallele\tdepth\tall\n";
while(my $line=<IN>){
	if($line =~ m/^#/){next;}	#Skip header lines
	chomp $line;
	my @data=split /\t/, $line;
	#Go through each genotype and output quality and genotype
	for(my $i=$startcol; $i < scalar @data; $i++){
		if($data[$i] eq "./."){next;}	#Skip missing data
		my @geno = split /:/, $data[$i];
		if($geno[$alleleID] eq "./."){next;}	#Skip missing genotypes (set by VCFtools)
		
		#Set genotype
		my $call;
		my @alleles = split /\//, $geno[$alleleID];
		if($alleles[0] ==$alleles[1]){
			$call="Hom";
		}else{
			$call="Het";
		}
		
		#print "$geno[$alleleID] to $call\n";
	
		#Get quality
		my $qual = $geno[$qualID];
		print OUT "$qual\t$call\t$geno[$depthID]\t$data[$i]\n";
	}
}
close IN;
close OUT;

exit;