# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $help = ''),
	'a=f' => \(my $alpha = 0.05),
	'D' => \(my $noDriveCount = ''),
	'O' => \(my $noOddsratio = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl VAMP_fisher.filter.pl VAMP_fisher.txt > VAMP_fisher.filter.txt

Options: -h       display this help message
         -a FLOAT alpha, p-value cutoff
         -D       do not consider drive count
         -O       do not consider odds ratio

EOF
}
my ($fisherFile) = @ARGV;
open(my $reader, $fisherFile);
my @columnList = ();
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ s/^#//) {
		@columnList = split(/\t/, $line);
		print '#', join("\t", @columnList), "\n";
	} else {
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		if($noDriveCount || $tokenHash{'driveCount'}) {
			if($tokenHash{'pvalue'} <= $alpha) {
				print join("\t", @tokenHash{@columnList}), "\n";
			} elsif($noOddsratio eq '' && ($tokenHash{'oddsratio'} eq "Inf" || $tokenHash{'oddsratio'} == 0)) {
				print join("\t", @tokenHash{@columnList}), "\n";
			}
		}
	}
}
close($reader);
