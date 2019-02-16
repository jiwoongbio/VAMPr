# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Scalar::Util qw(looks_like_number);
use List::Util qw(reduce any all none max min sum uniq);
use Statistics::R;

use Getopt::Long qw(:config no_ignore_case);

my @genotypeFileList = ();
GetOptions(
	'h' => \(my $help = ''),
	'g=s' => \@genotypeFileList,
	'c' => \(my $clusterOnly = ''),
	'o' => \(my $orthologyOnly = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl VAMP_fisher.pl [options] phenotype=VAMP.txt[,...] [...] > VAMP_fisher.txt

Options: -h       display this help message
         -g FILE  genotype file
         -c       cluster only
         -o       orthology only

EOF
}
my (@phenotypeFileList) = @ARGV;

my %genotypesHash = ();
foreach my $genotypeFile (@genotypeFileList) {
	open(my $reader, $genotypeFile);
	my @columnList = ();
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ s/^#//) {
			@columnList = split(/\t/, $line, -1);
		} elsif(@columnList) {
			my %tokenHash = ();
			@tokenHash{@columnList} = split(/\t/, $line, -1);
			$genotypesHash{$_} = 1 foreach(grep {defined} @tokenHash{'genotype', 'genotypes'});
		} else {
			my @tokenList = split(/\t/, $line, -1);
			$genotypesHash{$tokenList[0]} = 1;
		}
	}
	close($reader);
}
my @genotypesList = sort keys %genotypesHash;
my %genotypeHash = ();
$genotypeHash{$_} = 1 foreach(map {split(/,/, $_)} @genotypesList);

my $number = 0;
my %genotypeNumberHash = ();
my %genotypeClusterHash = ();
my %clusterNumberHash = ();
my %numberPhenotypeHash = ();
my %phenotypeHash = ();
my @phenotypeList = ();
foreach(map {[$_->[0], split(/,/, $_->[1])]} map {[split(/=/, $_, 2)]} @phenotypeFileList) {
	my ($phenotype, @fileList) = @$_;
	foreach my $file (@fileList) {
		$number += 1;
		open(my $reader, $file);
		while(my $line = <$reader>) {
			chomp($line);
			my ($query, $genotype) = split(/\t/, $line, -1);
			(my $cluster = $genotype) =~ s/\|.*$//;
			(my $orthology = $cluster) =~ s/\.[0-9]*$//;
			if($orthologyOnly) {
				$genotypeNumberHash{$orthology}->{$number} = 1;
			} elsif($clusterOnly) {
				$genotypeNumberHash{$cluster}->{$number} = 1;
			} else {
				if(@genotypeFileList) {
					$genotype = $cluster unless($genotypeHash{$genotype});
				}
				$genotypeNumberHash{$genotype}->{$number} = 1;
				$genotypeClusterHash{$genotype} = $cluster if($genotype ne $cluster);
				$clusterNumberHash{$cluster}->{$number} = 1;
			}
		}
		close($reader);
		$numberPhenotypeHash{$number} = $phenotype;
	}
	unless($phenotypeHash{$phenotype}) {
		$phenotypeHash{$phenotype} = 1;
		push(@phenotypeList, $phenotype) if($phenotype ne '');
	}
}

my %numbersGenotypeHash = ();
foreach my $genotype (sort keys %genotypeNumberHash) {
	my @numberList = sort {$a <=> $b} keys %{$genotypeNumberHash{$genotype}};
	my $numbers = join(',', @numberList);
	if(defined(my $cluster = $genotypeClusterHash{$genotype})) {
		my @numberList = sort {$a <=> $b} keys %{$clusterNumberHash{$cluster}};
		$numbers = join('|', $numbers, join(',', @numberList));
	}
	$numbersGenotypeHash{$numbers}->{$genotype} = 1 if(scalar(@numberList) < $number);
}

my %numbersColumnHash = ();
my %numberNumbersHash = ();
foreach my $numbers (keys %numbersGenotypeHash) {
	my @numberListList = map {[split(/,/, $_)]} split(/\|/, $numbers);
	my @numberList = @{$numberListList[0]};
	my %numberHash = map {$_ => 1} @numberList;
	$numberListList[0] = [1 .. $number];
	my @columnHashList = ();
	foreach(@numberListList) {
		my @numberList = @$_;
		my %phenotypeCountListHash = ();
		($_->[0], $_->[1]) = (0, 0) foreach(@phenotypeCountListHash{@phenotypeList});
		foreach my $number (@numberList) {
			my $phenotype = $numberPhenotypeHash{$number};
			if($numberHash{$number}) {
				$phenotypeCountListHash{$phenotype}->[0] += 1;
			} else {
				$phenotypeCountListHash{$phenotype}->[1] += 1;
			}
		}
		my @countList = map {($_->[0], $_->[1])} @phenotypeCountListHash{@phenotypeList};

		my %columnHash = ();
		@columnHash{map {("$_\_variant_count", "$_\_nonvariant_count")} @phenotypeList} = @countList;
		$columnHash{$_->[0]} = $_->[1] foreach(getFisherTestValueList(2, @countList));

		@columnHash{'phenotype', 'jaccardIndex'} = ($phenotypeList[0], $countList[0] / (sum(@countList) - $countList[3])) if($countList[0] * $countList[3] > $countList[2] * $countList[1]);
		@columnHash{'phenotype', 'jaccardIndex'} = ($phenotypeList[1], $countList[2] / (sum(@countList) - $countList[1])) if($countList[2] * $countList[1] > $countList[0] * $countList[3]);

		push(@columnHashList, \%columnHash);
	}

	my $pvalue = min(map {$columnHashList[$_]->{'pvalue'}} 0 .. $#columnHashList);
	my ($index) = grep {$columnHashList[$_]->{'pvalue'} eq $pvalue} 0 .. $#columnHashList;
	if($index == 0 && scalar(@numberListList) > 1) {
		my %genotypeHash = %{$numbersGenotypeHash{$numbers}};
		$numbers = join(',', @numberList);
		$numbersGenotypeHash{$numbers}->{$_} = 1 foreach(keys %genotypeHash);
	}
	unless(defined($numbersColumnHash{$numbers})) {
		$numbersColumnHash{$numbers} = $columnHashList[$index];
		$numbersColumnHash{$numbers}->{'mutationType'} = 'LOF' if($index > 0);
		if(defined(my $phenotype = $numbersColumnHash{$numbers}->{'phenotype'})) {
			foreach my $number (@{$numberListList[$index]}) {
				if($numberPhenotypeHash{$number} eq $phenotype) {
					$numberNumbersHash{$number}->{$numbers} = 1 if($numberHash{$number});
				} else {
					$numberNumbersHash{$number}->{$numbers} = 1 unless($numberHash{$number});
				}
			}
		}
	}
}
foreach my $number (keys %numberNumbersHash) {
	my @numbersList = keys %{$numberNumbersHash{$number}};
	my $pvalue = min(map {$numbersColumnHash{$_}->{'pvalue'}} @numbersList);
	foreach my $numbers (grep {$numbersColumnHash{$_}->{'pvalue'} eq $pvalue} @numbersList) {
		$numbersColumnHash{$numbers}->{'driveCount'} += 1;
	}
}

my @columnList = ('genotypes', 'pvalue', 'oddsratio', 'phenotype', 'jaccardIndex', 'driveCount', 'mutationType', map {("$_\_variant_count", "$_\_nonvariant_count")} @phenotypeList);
print '#', join("\t", @columnList), "\n";
foreach my $numbers (sort keys %numbersColumnHash) {
	$numbersColumnHash{$numbers}->{'genotypes'} = join(',', sort {compare([split(/\|/, $a)], [split(/\|/, $b)])} keys %{$numbersGenotypeHash{$numbers}});
	print join("\t", map {defined($_) ? $_ : ''} @{$numbersColumnHash{$numbers}}{@columnList}), "\n";
}

sub getFisherTestValueList {
	my ($nrow, @countList) = @_;
	my $R = Statistics::R->new();
	$R->run(sprintf('test.out <- fisher.test(matrix(c(%s), nrow = %d), simulate.p.value = TRUE)', join(',', @countList), $nrow));
	my $pvalue = $R->get('as.numeric(test.out$p.value)');
	my $oddsratio = $R->get('as.numeric(test.out$estimate["odds ratio"])');
	$R->stop();
	return (['pvalue', $pvalue], $oddsratio eq 'numeric(0)' ? () : ['oddsratio', $oddsratio]);
}

sub compare {
	my ($a, $b) = @_;
	my @a = @$a;
	my @b = @$b;
	if(scalar(@a) > 0 && scalar(@b) > 0) {
		$a = shift @a;
		$b = shift @b;
		if(looks_like_number($a) && looks_like_number($b)) {
			return $a <=> $b || compare(\@a, \@b);
		} else {
			return $a cmp $b || compare(\@a, \@b);
		}
	} elsif(scalar(@a) > 0) {
		return 1;
	} elsif(scalar(@b) > 0) {
		return -1;
	} else {
		return 0;
	}
}
