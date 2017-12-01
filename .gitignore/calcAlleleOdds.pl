#Considering the following sample results of a case-control study for a disease, write a Perl script 
#that study each position (assuming all positions are SNP positions)
#reporting the allele with the maximum odds ratio. 
#You need to calculate the odds ratio for each allele, then report the allele with the maximum
#odds ratio).
#In case, an allele with an odds ratio > 1.5 is found for a position, report that this position is
#associated with the disease.
#########################################################################################################
#!/usr/bin/perl
use strict;
use diagnostics;

#Command line arrguments for inputs case,control and output file alleles.tsv
my $inputCase = $ARGV[0];
my $inputControl = $ARGV[1];
my $oddsRatioOut = $ARGV[2];

#Hash Declaration
my %hashCase = ();
my %hashControl = ();
my %hashRatio = ();

#Open the output filehandle OUTFILE for writing
open(OUTFILE, ">", "$oddsRatioOut") || die "Cannot write into the output file $!";
#Printing the headers into the output file
print OUTFILE "Position\t", "Allele\t", "Odds_ratio\t", "Associated?\n";

#for and while loop to read the case file line by line, the length of each line being 46
for (my $i=0; $i <= 46; $i+=1) {
	open(INCASE, "$inputCase" ) || die "Major problem: cannot read the sam file $!";
	while ( my $inLine = <INCASE> ) {
#to memove end of line characters
		chomp($inLine);
#splitting into an array for easy access of induvidual bits
		my @inLineArray = split( "", $inLine );
#if loop to count the number of occurences of each base in each sequence of the case     
		if ($inLineArray[$i] eq 'A'){
			$hashCase{'A'} += 1;
		}
		elsif ($inLineArray[$i] eq 'C'){
			$hashCase{'C'} += 1;
		}
		elsif ($inLineArray[$i] eq 'T'){
			$hashCase{'T'} += 1;
		}
		elsif ($inLineArray[$i] eq 'G'){
			$hashCase{'G'} += 1;
		}

	}

	open(INCONTROL, "$inputControl" ) || die "Major problem: cannot read the sam file $!";
	while ( my $inLine = <INCONTROL> ) {
#to memove end of line characters
		chomp($inLine);
#splitting into an array for easy access of induvidual bits
		my @inLineArray = split( "", $inLine );
#if loop to count the number of occurences of each base in each sequence of the case
		if ($inLineArray[$i] eq 'A'){
			$hashControl{'A'} += 1;
		}
		elsif ($inLineArray[$i] eq 'C'){
			$hashControl{'C'} += 1;
		}
		elsif ($inLineArray[$i] eq 'T'){
			$hashControl{'T'} += 1;
		}
		elsif ($inLineArray[$i] eq 'G'){
			$hashControl{'G'} += 1;
		}

	}

	print OUTFILE "$i\t";
#calculating the required values for odds ratio
	foreach my $hashValue (keys %hashCase){
#presence of the alleles in the Cases is represented by the variable $a
		my $a = $hashCase{$hashValue};
#abscense of the alleles in the Cases is represented by the variable $c
		my $c = 8 - $hashCase{$hashValue};


		if (defined $hashControl{$hashValue}) {
#the number of controls with exposure is represented by the variable $b
			my $b = $hashControl{$hashValue};
#the number of controls without exposure is represented by the variable $d
			my $d = 8 - $hashControl{$hashValue};
# zeros cause problems with computation of the odds ratio so its standard error, 0.5 is added to all cells if any one cell is zero (a, b, c, d) 
			if ($a||$b||$c||$d == '0'){
				$a = ($a + 0.5);
				$c = ($c + 0.5);
				$b = ($b + 0.5);
				$d = ($d + 0.5);

				my $casesOddsExp = $a*$d;
				my $controlsOddsExp = $b*$c;

				my $oddsRatio = ($casesOddsExp)/($controlsOddsExp);
				$hashRatio{$hashValue} = $oddsRatio;

			}
			else {
				my $casesOddsExp = $a*$d;
				my $controlsOddsExp = $b*$c;

				my $oddsRatio = ($casesOddsExp)/($controlsOddsExp);
				$hashRatio{$hashValue} = $oddsRatio;
			}
		}

		elsif (not defined $hashControl{$hashValue}) {

			my $b = 0;
			my $d = 8 - 0;

			if ($a||$b||$c||$d == '0'){
				$a = ($a + 0.5);
				$c = ($c + 0.5);
				$b = ($b + 0.5);
				$d = ($d + 0.5);

				my $casesOddsExp = $a*$d;
				my $controlsOddsExp = $b*$c;

				my $oddsRatio = ($casesOddsExp)/($controlsOddsExp);
				$hashRatio{$hashValue} = $oddsRatio;

			}
			else {
				my $casesOddsExp = $a*$d;
				my $controlsOddsExp = $b*$c;

				my $oddsRatio = ($casesOddsExp)/($controlsOddsExp);
				$hashRatio{$hashValue} = $oddsRatio;
			}
		} 

	}

	my @sortedhashRatio = reverse sort keys %hashRatio;

	if (exists $sortedhashRatio[1]){
		if ($hashRatio{$sortedhashRatio[0]} == $hashRatio{$sortedhashRatio[1]}) {
			print OUTFILE "$sortedhashRatio[0],$sortedhashRatio[1]\t", "$hashRatio{$sortedhashRatio[0]}\t";
		}
		elsif ($hashRatio{$sortedhashRatio[0]} > $hashRatio{$sortedhashRatio[1]}){
			print OUTFILE "$sortedhashRatio[0]\t", "$hashRatio{$sortedhashRatio[0]}\t";
		}
	}
	else {
		print OUTFILE "$sortedhashRatio[0]\t", "$hashRatio{$sortedhashRatio[0]}\t";
	}



	if ($hashRatio{$sortedhashRatio[0]} > '1.5'){
		print OUTFILE "Y\n";
	}
	else{
		print OUTFILE "N\n";
	}


	%hashCase = ();
	%hashControl =();
	%hashRatio = ();

}

close(INCASE);
close(INCONTROL);

