#!/usr/bin/env perl

#usage: $PSHOME/scripts/second_allele_calculations.pl $race counts1.R0k6 $ids 1 $PSHOME $outDir


use Bio::DB::Sam;
use List::MoreUtils 'true';
use List::MoreUtils 'indexes';
use List::MoreUtils qw/ uniq /;
use List::Util qw/max/;
#use FindBin qw($Bin);
#use lib "$Bin/../lib";
use PolysolverUtils qw(get_hla_frequencies population_prior_likelihood_component);

$race = $ARGV[0];
$winners1File = $ARGV[1];
$idsFile = $ARGV[2];
$freqFile = $ARGV[3];
$a1_dir = $ARGV[4];
$outDir = $ARGV[5];

my $includeFreq = 1;

$scale = exp(23);
#$freqFile = $PSHOME."/"."data/HLA_FREQ.txt";
$offset = 0;
$scoreCol = 3;

if ( ! -d $outDir ) {
	mkdir($outDir);
}
# get hla frequencies

%freqHash = %{get_hla_frequencies($freqFile,"f4")};

# get population frequency component

# FIXME: sort causes a series of broken pipe errors here
# $alleleA1 = `sort -k3,2rn $winners1File | grep hla_a | head -1 | cut -f1`;
$alleleA1 = `sort -k2,2nr $winners1File | grep hla_a | head -n1 | cut -f1`;
chomp($alleleA1);
$alleleB1 = `sort -k2,2rn $winners1File | grep hla_b | head -n1 | cut -f1`;
chomp($alleleB1);
$alleleC1 = `sort -k2,2rn $winners1File | grep hla_c | head -n1 | cut -f1`;
chomp($alleleC1);

print "winners1\t$alleleA1\t$alleleB1\t$alleleC1\n";

# get read likelihoods for winners in hash

%winnerA = ();
%winnerB = ();
%winnerC = ();

%winnerA = %{get_winner1_allele_reads($alleleA1,$includeFreq,$scoreCol,$offset,$a1_dir)};
%winnerB = %{get_winner1_allele_reads($alleleB1,$includeFreq,$scoreCol,$offset,$a1_dir)};
%winnerC = %{get_winner1_allele_reads($alleleC1,$includeFreq,$scoreCol,$offset,$a1_dir)};

# get modified likelihood for given allele (allele==allele)

open IDSFILE, $idsFile;

$count = 0;

while($allele = <IDSFILE>){

	chomp($allele);

	$count++;

	#print "$count\t$allele\n";
	$freqComponent = 0;
	$freqComponentLog = 0;

	if($race eq "Caucasian" || $race eq "Black" || $race eq "Asian"){
        	$freqComponent = population_prior_likelihood_component(\%freqHash,$allele,$race,$scale);
		if($freqComponent > 0){
     		   	$freqComponentLog = log($freqComponent);
		}
	}

	if($includeFreq == 1){
		$inFile = $a1_dir."/".$allele.".lik1";
		$outFile = $outDir."/".$allele.".lik2";
	} else {
		$inFile = $a1_dir."/".$allele.".nofreq.lik1";
		$outFile = $outDir."/".$allele.".nofreq.lik2";
	}

	open INFILE, $inFile;
	
	open OUTFILE, ">$outFile";

	$likLogIscoreTotal = 0;
	$reads = 0;
	while($line = <INFILE>){
		chomp($line);

		if($line =~ /^lik1/){
			last;
		}

		$reads = $reads + 1;

		@f = split(/\t/,$line);

		my $pair = $f[0];

		my $l1;
		my $l2 = $f[$scoreCol-1] + $offset;

		$iScoreLog = $f[$scoreCol];

		if($allele =~ /hla_a/){
			if(exists $winnerA{$pair}) {
				$l1 = $winnerA{$pair};
			}else {
				$l1 = 0;
			}
		} elsif($allele =~ /hla_b/){
			if(exists $winnerB{$pair}) {
				$l1 = $winnerB{$pair};
			}else {
				$l1 = 0;
			}
		} elsif($allele =~ /hla_c/){
			if(exists $winnerC{$pair}) {
				$l1 = $winnerC{$pair};
			}else {
				$l1 = 0;
			}
    } else {
			die "allele = $allele is not valid in parse_sam_lik2_fast.pl pair = $pair\n";
		}
		
		if(abs($l1) > 0 & abs($l2) > 0){
			$factor = $l2/($l1+$l2);
			$l2 = $factor*$l2 + $iScoreLog;
			#$l2 = 0.5;
		}		

		#print "$pair\t$factor\t$l2\n";
		# Note: I dont need to output per-reads score for 2nd allele
		# only needs the final score for typing
		#print OUTFILE "$pair\t$l2\n";

		$likLogIscoreTotal = $likLogIscoreTotal + $l2;
	}

	if($includeFreq == 1){
    $likLogIscoreTotalFreq = $likLogIscoreTotal + $freqComponentLog;
	} else {
    $likLogIscoreTotalFreq = $likLogIscoreTotal;
	}

	print OUTFILE "lik2\t$likLogIscoreTotalFreq\t$likLogIscoreTotal\t$freqComponentLog\n";

	close INFILE;

	close OUTFILE;

}

close IDSFILE;


########

sub get_winner1_allele_reads { # start sub get_winner1_allele_reads

	my ($name,$includeFreq,$scoreCol,$offset,$outDir) = @_;

	my ($line,%hash,@f,$inFile);

	if($includeFreq == 1){
        	$inFile = $outDir."/".$name.".lik1";
	} else{
        	$inFile = $outDir."/".$name.".nofreq.lik1";
	}
	if ( ! -e $inFile) {
		die "Cannot find the $inFile\n";
	}

	open INFILE, $inFile;
	while($line = <INFILE>){
 		chomp($line);
 		if($line=~/^lik/){
        		last;
 		}
 		@f = split(/\t/,$line);
		$hash{$f[0]} = $f[$scoreCol - 1];
	}
	close INFILE;

	return \%hash;

} # end sub get_winner1_allele_reads
