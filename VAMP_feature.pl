# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Cwd 'abs_path';
use Getopt::Long qw(:config no_ignore_case);

(my $vampPath = abs_path($0)) =~ s/\/[^\/]*$//;
my $dataPath = "$vampPath/VAMP_data";

GetOptions('h' => \(my $help = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl VAMP_feature.pl [options] VAMP.txt > VAMP_feature.txt

Options: -h       display this help message

EOF
}
my ($file) = @ARGV;

my %clusterStartEndTokenListListHash = ();
{
	open(my $reader, $file);
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line, -1);
		my ($query, $genotype) = @tokenList;
		my ($cluster, $clusterPosition, $clusterAA, $variantAA) = split(/\|/, $genotype);
		if(defined($clusterPosition)) {
			my ($genotypeStart, $genotypeEnd) = ($clusterPosition, $clusterPosition + length($clusterAA) - 1);
			push(@{$clusterStartEndTokenListListHash{$cluster}}, [$genotypeStart, $genotypeEnd, @tokenList]);
		}
	}
	close($reader);
}

{
	open(my $reader, "$dataPath/cluster.feature.txt");
	while(my $line = <$reader>) {
		chomp($line);
		my ($cluster, $start, $end, $clusterSequence, $original, $variation, $type, $description) = split(/\t/, $line, -1);
		if(defined(my $startEndTokenListList = $clusterStartEndTokenListListHash{$cluster})) {
			foreach(@$startEndTokenListList) {
				my ($genotypeStart, $genotypeEnd, @tokenList) = @$_;
				if(($start eq '' && $end eq '') || ($start <= $genotypeEnd && $genotypeStart <= $end)) {
					print join("\t", @tokenList, $cluster, $start, $end, $clusterSequence, $original, $variation, $type, $description), "\n";
				}
			}
		}
	}
	close($reader);
}
