# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(reduce any all none max min sum uniq);
use Statistics::R;

use Getopt::Long qw(:config no_ignore_case);

my @phenotypeList = ();
GetOptions(
	'h' => \(my $help = ''),
	'c' => \(my $clusterOnly = ''),
	'o' => \(my $orthologyOnly = ''),
	'p=s' => \@phenotypeList,
	'i=s' => \(my $imageFile = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl VAMP_predict.pl [options] VAMP_model.RData [sample=]VAMP.txt [...]

Options: -h       display this help message
         -c       cluster only
         -o       orthology only
         -p STR   phenotype
         -i FILE  image file

EOF
}
my ($modelFile, @sampleFileList) = @ARGV;
if(-r $modelFile) {
	my $R = Statistics::R->new();
	$R->run(sprintf('load("%s")', $modelFile));
	$R->run('library(caret)');
	$R->run('library(xgboost)');
	my @genotypesList = grep {$_ ne '.outcome'} @{$R->get('colnames(model$trainingData)')};
	my %genotypeHash = ();
	$genotypeHash{$_} = 1 foreach(map {split(/,/, $_)} @genotypesList);

	my @sampleList = ();
	my %sampleGenotypeHash = ();
	foreach(map {[$_->[0], $_->[-1]]} map {[split(/=/, $_, 2)]} @sampleFileList) {
		my ($sample, $file) = @$_;
		push(@sampleList, $sample);
		open(my $reader, $file);
		while(my $line = <$reader>) {
			chomp($line);
			my ($query, $genotype) = split(/\t/, $line, -1);
			(my $cluster = $genotype) =~ s/\|.*$//;
			(my $orthology = $cluster) =~ s/\.[0-9]*$//;
			if($orthologyOnly) {
				$sampleGenotypeHash{$sample}->{$orthology} = 1;
			} elsif($clusterOnly) {
				$sampleGenotypeHash{$sample}->{$cluster} = 1;
			} else {
				$genotype = $cluster unless($genotypeHash{$genotype});
				$sampleGenotypeHash{$sample}->{$genotype} = 1;
			}
		}
		close($reader);
	}

	$R->run('x <- data.frame()');
	foreach my $index (0 .. $#sampleList) {
		my $sample = $sampleList[$index];
		my %genotypeHash = %{$sampleGenotypeHash{$sample}};
		$R->run(sprintf('x <- rbind(x, matrix(c(%s), nrow = 1))', join(',', map {(all {$genotypeHash{$_}} split(/,/, $_)) ? 1 : 0} @genotypesList)));
	}
	$R->set('genotypes', \@genotypesList);
	$R->run('colnames(x) <- genotypes');
	$R->run('x <- data.matrix(x)');
	$R->run('y <- predict(model, x)');
	$R->run('probs <- predict(model, x, type = "prob")');
	foreach my $index (0 .. $#sampleList) {
		my $sample = $sampleList[$index];
		my $phenotype = $R->get(sprintf('as.character(y[%d])', $index + 1));
		my @probabilityList = map {$R->get(sprintf('as.numeric(probs[%d, "%s"])', $index + 1, $_))} @phenotypeList;
		print join("\t", $sample, $phenotype, @probabilityList), "\n";
	}
	$R->run(sprintf('save.image(file = "%s")', $imageFile)) if($imageFile ne '');
	$R->stop();
}
