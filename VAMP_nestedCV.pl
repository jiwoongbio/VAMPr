# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Scalar::Util qw(looks_like_number);
use List::Util qw(reduce any all none max min sum uniq);
use Statistics::R;

use Getopt::Long qw(:config no_ignore_case);

my @genotypeFileList = ();
GetOptions('h' => \(my $help = ''),
	'g=s' => \@genotypeFileList,
	'c' => \(my $clusterOnly = ''),
	'o' => \(my $orthologyOnly = ''),
	's=i' => \(my $seed = 1),
	'f=i' => \(my $fold = 5),
	'F=i' => \(my $outerFold = 10),
	'G=s' => \(my $grid = ''),
	'i=s' => \(my $imageFile = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl VAMP_nestedCV.pl [options] VAMP_nestedCV.accuracy.txt phenotype=VAMP.txt,... [...]

Options: -h       display this help message
         -g FILE  genotype file
         -c       cluster only
         -o       orthology only
         -s INT   seed [$seed]
         -f INT   fold [$fold]
         -F INT   outer fold [$outerFold]
         -G STR   grid
         -i FILE  image file

EOF
}
my ($accuracyFile, @phenotypeFileList) = @ARGV;

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
my %numberPhenotypeHash = ();
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
			}
		}
		close($reader);
		$numberPhenotypeHash{$number} = $phenotype;
	}
}

my %numberGenotypesHash = ();
if(@genotypeFileList) {
	foreach my $genotypes (@genotypesList) {
		my @genotypeList = split(/,/, $genotypes);
		my %numberGenotypeHash = ();
		foreach my $genotype (@genotypeList) {
			$numberGenotypeHash{$_}->{$genotype} = 1 foreach(keys %{$genotypeNumberHash{$genotype}});
		}
		foreach my $number (keys %numberGenotypeHash) {
			my %genotypeHash = %{$numberGenotypeHash{$number}};
			$numberGenotypesHash{$number}->{$genotypes} = 1 if(all {$genotypeHash{$_}} @genotypeList);
		}
	}
} else {
	my %numbersGenotypeHash = ();
	foreach my $genotype (sort keys %genotypeNumberHash) {
		my @numberList = sort {$a <=> $b} keys %{$genotypeNumberHash{$genotype}};
		my $numbers = join(',', @numberList);
		$numbersGenotypeHash{$numbers}->{$genotype} = 1 if(scalar(@numberList) < $number);
	}
	foreach my $numbers (sort keys %numbersGenotypeHash) {
		my $genotypes = join(',', sort {compare([split(/\|/, $a)], [split(/\|/, $b)])} keys %{$numbersGenotypeHash{$numbers}});
		push(@genotypesList, $genotypes);
		foreach my $number (split(/,/, $numbers)) {
			$numberGenotypesHash{$number}->{$genotypes} = 1;
		}
	}
}

{
	my $R = Statistics::R->new();
	$R->run('x <- data.frame()');
	$R->run('y <- c()');
	foreach my $number (1 .. $number) {
		my %genotypesHash = defined($_ = $numberGenotypesHash{$number}) ? %$_ : ();
		foreach my $index (0 .. $#genotypesList) {
			$R->set(sprintf('x[%d, %d]', $number, $index + 1), $genotypesHash{$genotypesList[$index]} ? 1 : 0);
		}
		$R->set("y[$number]", $numberPhenotypeHash{$number});
	}
	foreach my $index (0 .. $#genotypesList) {
		$R->run('genotypes <- c()');
		my @genotypeList = split(/,/, $genotypesList[$index]);
		foreach my $index (0 .. $#genotypeList) {
			$R->set(sprintf('genotypes[%d]', $index + 1), $genotypeList[$index]);
		}
		$R->run(sprintf('colnames(x)[%d] <- paste(genotypes, collapse = ",")', $index + 1));
	}
	$R->run('x <- data.matrix(x)');
	$R->run('y <- factor(y)');
	$R->run('library(caret)');
	$R->run('library(xgboost)');
	$R->run("set.seed($seed)");
	$R->run("folds <- createFolds(y, k = $outerFold, list = FALSE)");
	$R->run("models <- list()");
	$R->run(sprintf('for(i in unique(folds)) {if(length(unique(y[folds != i])) == length(levels(y))) {models[[i]] <- train(x[folds != i, ], y[folds != i], %s)}}', join(', ',
		'method = "xgbTree"',
		sprintf('trControl = trainControl(method = "repeatedcv", number = %d, repeats = 1, classProbs = TRUE, allowParallel= TRUE)', $fold),
		($grid eq '' ? () : "tuneGrid = expand.grid($grid)"),
		'metric = "Accuracy"',
	)));
	$R->run(sprintf('write.table(t(sapply(unique(folds), function(i) {c(%s)})), file = "%s", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)', join(', ',
		'ifelse(is.null(models[[i]]), NA, mean(predict(models[[i]], x[folds == i, ]) == y[folds == i]))',
		'ifelse(is.null(models[[i]]), NA, max(models[[i]]$results$Accuracy))',
	), $accuracyFile));
	$R->run(sprintf('save.image(file = "%s")', $imageFile)) if($imageFile ne '');
	$R->stop();
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
