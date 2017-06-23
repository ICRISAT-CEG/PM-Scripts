#! /usr/bin/perl -w

#Prupose: Reformat TASSEL's hapmap output for processing with MSTmap
#Arguments: (0) Input file name (hapmap file), (1) Output file name, (2) Parent1 names (comma-separated), (3) Parent2 names (comma-separated),
#	(4) Population type (RILx, DH, etc; DH is for whenever you don't expect hets), (5) p-value cutoff (<1 will separate into different groups; >1 will assume all 1 group)
#Note: currenlty hardcoded for RILs; minor modifications for others. Also still has some print statements to aid debugging

use strict;

if(scalar(@ARGV)!=8){print "Must supply arguments of input file, output file, parent1, parent2, poptype, pval, min LG inclusion distance, and min LG size cutoff\n"; exit;}
my $infile=$ARGV[0];
my $outfile=$ARGV[1];
my $parent1=$ARGV[2];
my $parent2=$ARGV[3];
my $poptype=$ARGV[4];
my $pval=$ARGV[5];
my $nomap_dist=$ARGV[6];
my $nomap_size=$ARGV[7];

print "Converting $infile to MST format\n";

#columns in hapmap file with data
my $snpID=0;
my $allelesID=1;
my $startID=11;

#Record parental names
my %p1;
my %p2;
foreach my $p (split /,/, $parent1){
	$p1{$p}=1;
}
foreach my $p (split /,/, $parent2){
	$p2{$p}=1;
}

open(INFILE, "$infile") or die "Cannot open file $infile.\n";
my $header=<INFILE>;
chomp $header;
#search through header for position of parents
my @headerdata = split/\t/, $header;
my @par1;
my @par2;
my $needToMerge=0;
my %toskip=();	#hash of locations where parents are
for(my $n=$startID; $n<scalar @headerdata; $n++){
	#$headerdata[$n] =~ /(\w+)/;	#capture first wordlike part of name
	#my $name=$1;
	my $name = $headerdata[$n];
	#print "Name is $name\n";
	if (exists $p1{$name}){push @par1, $n; $toskip{$n}=1;}
	elsif (exists $p2{$name}){push @par2, $n; $toskip{$n}=1;}
}
if(scalar @par1>1 || scalar @par2>1){$needToMerge=1;}	#flag for if need to merge multiple parents.


#Go through file and convert, one line at a time
my $numloci=0;
my @loci;
while (my $line=<INFILE>){
	chomp $line;
	my @data = split /\t/, $line;
	my @alleles = split /\//, $data[$allelesID];
	my $locus="$data[$snpID]";
	#assign parental genotype
	my $p1 = $data[$par1[0]];
	my $p2 = $data[$par2[0]];
	#print "Data start as $p1 and $p2\n";
	if ($needToMerge==1){
		my $p1alleles = $data[$par1[0]];		#INCORRECT JOIN; NEED TO FIX
		for(my $n=1; $n< scalar @par1; $n++){$p1alleles .= $data[$par1[$n]];}
		$p1 = merge($p1alleles);
		my $p2alleles = $data[$par2[0]];
		for(my $n=1; $n< scalar @par2; $n++){$p2alleles .= $data[$par2[$n]];}
		$p2 = merge($p2alleles);
		print "\tAllele calls $p1alleles $p2alleles\n";
	}
	#print "\tData end as as $p1 and $p2\n";
	if($p1 eq $p2){print "\tSkipping because the same $p1 and $p2\n"; next;}	#if parents supposedly have same genotype, skip
	#really crude imputation; assume a missing parent has the genotype the other doesn't	SEE IF THIS ACTUALLY HELPS OR NOT
	if($p1 eq "N"){
		if($p2 eq $alleles[0]){$p1=$alleles[1];}
		else {$p1=$alleles[0];}
	}
	if($p2 eq "N"){
		if($p1 eq $alleles[0]){$p2=$alleles[1];}
		else {$p2=$alleles[0];}
	}
	
	#go through and assign genotypes
	for(my $n=$startID; $n<scalar @data; $n++){
		#if(exists $toskip{$n}){next;}	#skip genotype if is a parental line
		my $call;
		if($data[$n] eq "N"){$call = "U";}	#missing data
		elsif($data[$n] =~ m/[RYSWKMH]/) {	#heterozygote - set to missing if pop shouldn't have any, or to het if is a RIL
			if($poptype eq "DH") {$call="U";}
			else{$call="X";}
		}	
		elsif($data[$n] eq $p1){$call = "A";}	#parent1
		elsif($data[$n] eq $p2){$call = "B";}	#parent2
		else {$call="U";}	#default to unknown
		$locus .= "\t$call";
		#if($n - $startID % 5==0){$locus .= " ";}	#add spaces every five calls
		#print "\tAllele $data[$n] turned to $call for parents $p1 and $p2\n";
		#print "\tProcessing $headerdata[$n] to $call\n";
	}
	$locus.="\n";
	push @loci, $locus;
	$numloci++;
}
close INFILE;

open (OUTFILE, ">$outfile") or die "Cannot open file $outfile.\n";
#print header
$infile =~ /(\w+)\./;	#capture file name with regex and assign to population name
my $popname =$1;
my $numlines = scalar(@headerdata) - $startID;# - scalar(@par1) - scalar(@par2);
print OUTFILE "population_type $poptype\n";
print OUTFILE "population_name $popname\n";
print OUTFILE "distance_function kosambi\n";
print OUTFILE "cut_off_p_value $pval\n";
print OUTFILE "no_map_dist $nomap_dist\n";				#groups more than this cM from others are considered bad if have <[below] SNPs in them
print OUTFILE "no_map_size $nomap_size\n";				#groups less than this are considered bad if >[above] cM from others
print OUTFILE "missing_threshold 0.75\n";		#fraction mission data that removes marker (1.0= only if all missing)
print OUTFILE "estimation_before_clustering no\n";
print OUTFILE "detect_bad_data yes\n";
print OUTFILE "objective_function COUNT\n";	
print OUTFILE "number_of_loci $numloci\n";
print OUTFILE "number_of_individual $numlines\n\n";
print OUTFILE "locus_name";
for(my $n=$startID; $n<scalar @headerdata; $n++){
	#$headerdata[$n] =~ /(\w+)/;	#capture first wordlike part of name
	#my $name=$1;
	my $name = $headerdata[$n];
	#unless (exists $p1{$name} || exists $p2{$name}){
		print OUTFILE "\t$name";
		#}
}
print OUTFILE "\n";
print OUTFILE @loci;


#subroutine to merge multiple parental lines into one
sub merge{
	my $arg=shift;
	my @alleles = split //, $arg;
	for(my $n=1; $n<scalar @alleles; $n++){
		if($alleles[$n] eq $alleles[0]){next;}	#if both are equal, go to next allele
		elsif($alleles[$n] eq "N"){next;}			#if missing data, go to next allele
		elsif($alleles[0] eq "N"){$alleles[0] = $alleles[$n]; next;}	#if first allele is N and second isn't, reassign first as second
		else {$alleles[0] = "h";}		#if get past all else (by them not matching and neither being missing), then assign as heterozygote
	}
	return $alleles[0];
}

exit;