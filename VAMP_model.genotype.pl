# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);
use Statistics::R;

GetOptions('h' => \(my $help = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl VAMP_model.genotype.pl [options] VAMP_model.RData > VAMP_model.genotype.txt

Options: -h       display this help message

EOF
}
my ($modelFile) = @ARGV;

if(-r $modelFile) {
	my $R = Statistics::R->new();
	$R->run(sprintf('load("%s")', $modelFile));
	$R->run('library(caret)');
	$R->run('library(xgboost)');
	$R->run('importance <- varImp(model)$importance');
	my @variablesList = ref($_ = $R->get('rownames(importance)')) ? @$_ : ($_);
	my @importanceList = ref($_ = $R->get('importance$Overall')) ? @$_ : ($_);
	foreach my $index (0 .. $#variablesList) {
		print join("\t", $variablesList[$index], $importanceList[$index]), "\n";
	}
	$R->stop();
}
