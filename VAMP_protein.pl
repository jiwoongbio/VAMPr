# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Cwd 'abs_path';
use Getopt::Long qw(:config no_ignore_case);
use List::Util qw(max sum);
use Bio::DB::Fasta;

my %aaTypeHash = (
	'R' => 'b', 'H' => 'b', 'K' => 'b', # basic (b)
	'D' => 'a', 'E' => 'a', # acidic (a)
	'S' => 'p', 'T' => 'p', 'N' => 'p', 'Q' => 'p', # polar (p)
	'A' => 'h', 'V' => 'h', 'I' => 'h', 'L' => 'h', 'M' => 'h', 'F' => 'h', 'Y' => 'h', 'W' => 'h', # hydrophobic (h)
);

my @samMandatoryFieldList = ('qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual');

(my $vampPath = abs_path($0)) =~ s/\/[^\/]*$//;
my $dataPath = "$vampPath/VAMP_data";

GetOptions('h' => \(my $help = ''),
	't=s' => \(my $temporaryDirectory = defined($ENV{'TMPDIR'}) ? $ENV{'TMPDIR'} : '/tmp'),
	'p=i' => \(my $threads = 1),
	'e=f' => \(my $evalue = 10),
	'c=f' => \(my $minimumCoverage = 0.8),
	's=s' => \(my $samFile = ''),
	'a=s' => \(my $alignmentFile = ''),
	'A' => \(my $allVariants = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl VAMP_protein.pl [options] protein.fasta > VAMP.txt

Options: -h       display this help message
         -t DIR   directory for temporary files [\$TMPDIR or /tmp]
         -p INT   number of threads [$threads]
         -e FLOAT maximum e-value to report alignments [$evalue]
         -c FLOAT minimum coverage [$minimumCoverage]
         -s FILE  output SAM file
         -a FILE  output alignment file
         -A       all variants

EOF
}
my ($inputFastaFile) = @ARGV;
my %proteinSequenceHash = ();
{
	my $protein = '';
	open(my $reader, $inputFastaFile);
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^>(\S*)/) {
			$protein = $1;
		} else {
			$proteinSequenceHash{$protein} .= uc($line);
		}
	}
	close($reader);
}
chomp(my $hostname = `hostname`);
my $temporaryPrefix = "$temporaryDirectory/VAMP.$hostname.$$";
{
#	system("diamond blastp --threads $threads --db $dataPath/protein.dmnd --query $inputFastaFile --out $temporaryPrefix.sam --outfmt 101 --max-target-seqs 0 --evalue $evalue --unal 0 --tmpdir $temporaryDirectory --quiet");
	system("diamond blastp --threads $threads --db $dataPath/protein.dmnd --query $inputFastaFile --out $temporaryPrefix.sam --outfmt 101 --top 0 --evalue $evalue --unal 0 --tmpdir $temporaryDirectory --quiet");
}
{
	open(my $reader, "$temporaryPrefix.sam");
	open(my $writer, "| sort -t'\t' -k4,4 -k1,1g -k2,2gr -k3,3gr | cut -f4- > $temporaryPrefix.filtered.sam");
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^@/);
		my %tokenHash = ();
		(@tokenHash{@samMandatoryFieldList}, my @tagTypeValueList) = split(/\t/, $line);
		$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
		next if($tokenHash{'flag'} & 4);
		next if($tokenHash{'flag'} & 16);
		if(($tokenHash{'ZC:f'} = sum(0, $tokenHash{'MD:Z'} =~ /([0-9]+)/g) / $tokenHash{'ZL:i'}) >= $minimumCoverage) {
			push(@tagTypeValueList, "ZC:f:$tokenHash{'ZC:f'}");
			print $writer join("\t", @tokenHash{'ZE:f', 'AS:i', 'ZC:f'}, @tokenHash{@samMandatoryFieldList}, @tagTypeValueList), "\n";
		}
	}
	close($reader);
	close($writer);
}
{
	open(my $reader, "$temporaryPrefix.filtered.sam");
	open(my $writer, "> $temporaryPrefix.bestfit.sam");
	my %topTokenHash = ();
	$topTokenHash{'qname'} = '';
	while(my $line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		(@tokenHash{@samMandatoryFieldList}, my @tagTypeValueList) = split(/\t/, $line);
		$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
		%topTokenHash = %tokenHash if($tokenHash{'qname'} ne $topTokenHash{'qname'});
		if(scalar(grep {$tokenHash{$_} != $topTokenHash{$_}} 'ZE:f', 'AS:i', 'ZC:f') == 0) {
			@tagTypeValueList = map {join(':', $_, $tokenHash{$_})} map {join(':', @$_[0, 1])} map {[split(/:/, $_, 3)]} @tagTypeValueList;
			print $writer join("\t", @tokenHash{@samMandatoryFieldList}, @tagTypeValueList), "\n";
		}
	}
	close($reader);
	close($writer);
}
my @alignmentTokenListList = ();
{
	{
		my $db = Bio::DB::Fasta->new("$dataPath/cluster.fasta");
		sub getClusterSequence {
			my ($cluster, $start, $end) = @_;
			return $db->seq($cluster, $start, $end);
		}
	}
	open(my $reader, "$temporaryPrefix.bestfit.sam");
	open(my $writer, "| sort -u");
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^@/);
		my %tokenHash = ();
		(@tokenHash{@samMandatoryFieldList}, my @tagTypeValueList) = split(/\t/, $line);
		$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
		next if(scalar(my ($protein, $cluster, $clusterStart, $clusterEnd, $clusterCigar) = split(/\|/, $tokenHash{'rname'})) < 5);
		my @proteinIndexList = map {defined($_) ? $_ - 1 : $_} getPositionList(@tokenHash{'pos', 'cigar'});
		my @clusterIndexList = map {defined($_) ? $_ - 1 : $_} getPositionList(1, $clusterCigar);
		my @clusterVariantIndexList = grep {defined($_->[0]) && defined($_->[0] = $clusterIndexList[$_->[0]])} map {[$proteinIndexList[$_], $_ + ($tokenHash{'ZS:i'} - 1)]} 0 .. $#proteinIndexList;
		my $clusterSequence = getClusterSequence($cluster, $clusterStart, $clusterEnd);
		my $variantSequence = $proteinSequenceHash{$tokenHash{'qname'}};
		my @clusterAAList = split(//, $clusterSequence);
		my @variantAAList = split(//, $variantSequence);
		my @clusterAlignmentAAList = ();
		my @variantAlignmentAAList = ();
		my $alignmentLength = 0;
		{
			my $clusterNextIndex = 0;
			my $variantNextIndex = 0;
			foreach(@clusterVariantIndexList) {
				my ($clusterIndex, $variantIndex) = @$_;
				if($clusterIndex > $clusterNextIndex) {
					my $length = scalar(my @indexList = $clusterNextIndex .. $clusterIndex - 1);
					push(@clusterAlignmentAAList, @clusterAAList[@indexList]);
					push(@variantAlignmentAAList, ('-') x $length);
					$alignmentLength += $length;
				}
				if($variantIndex > $variantNextIndex) {
					my $length = scalar(my @indexList = $variantNextIndex .. $variantIndex - 1);
					push(@clusterAlignmentAAList, ('-') x $length);
					push(@variantAlignmentAAList, @variantAAList[@indexList]);
					$alignmentLength += $length;
				}
				push(@clusterAlignmentAAList, $clusterAAList[$clusterIndex]);
				push(@variantAlignmentAAList, $variantAAList[$variantIndex]);
				$clusterNextIndex = $clusterIndex + 1;
				$variantNextIndex = $variantIndex + 1;
				$alignmentLength += 1;
			}
			if((my $clusterIndex = scalar(@clusterAAList)) > $clusterNextIndex) {
				my $length = scalar(my @indexList = $clusterNextIndex .. $clusterIndex - 1);
				push(@clusterAlignmentAAList, @clusterAAList[@indexList]);
				push(@variantAlignmentAAList, ('-') x $length);
				$alignmentLength += $length;
			}
			if((my $variantIndex = scalar(@variantAAList)) > $variantNextIndex) {
				my $length = scalar(my @indexList = $variantNextIndex .. $variantIndex - 1);
				push(@clusterAlignmentAAList, ('-') x $length);
				push(@variantAlignmentAAList, @variantAAList[@indexList]);
				$alignmentLength += $length;
			}
		}
		if($alignmentFile ne '') {
			push(@alignmentTokenListList, [$tokenHash{'qname'}, $cluster, join('', @variantAlignmentAAList), join('', @clusterAlignmentAAList), $clusterStart]);
		}
		my @clusterPositionList = ();
		my %clusterPositionIndexHash = ();
		{
			my $clusterPosition = $clusterStart;
			foreach my $index (0 .. $alignmentLength - 1) {
				$clusterPositionList[$index] = $clusterPosition;
				if($clusterAlignmentAAList[$index] ne '-') {
					$clusterPositionIndexHash{$clusterPosition} = $index;
					$clusterPosition += 1;
				}
			}
		}
		my @variantPositionList = ();
		my %variantPositionIndexHash = ();
		{
			my $variantPosition = 1;
			foreach my $index (0 .. $alignmentLength - 1) {
				$variantPositionList[$index] = $variantPosition;
				if($variantAlignmentAAList[$index] ne '-') {
					$variantPositionIndexHash{$variantPosition} = $index;
					$variantPosition += 1;
				}
			}
		}
		my @genotypeList = ();
		foreach my $index (0 .. $alignmentLength - 1) {
			if(isMatched($clusterAlignmentAAList[$index], $variantAlignmentAAList[$index])) {
				$clusterAlignmentAAList[$index] = '=';
				$variantAlignmentAAList[$index] = '=';
			}
		}
		my $clusterAlignment = join('', @clusterAlignmentAAList);
		while($clusterAlignment =~ /([^=]+)/g) {
			my ($startIndex, $endIndex) = (pos($clusterAlignment) - length(my $clusterAA = $1), pos($clusterAlignment) - 1);
			my $variantAA = join('', @variantAlignmentAAList[$startIndex .. $endIndex]);
			my $clusterPosition = $clusterPositionList[$startIndex];
			my $variantPosition = $variantPositionList[$startIndex];
			$clusterAA =~ s/-//g;
			$variantAA =~ s/-//g;
			my $terminated = ($endIndex == $alignmentLength - 1 && ($clusterAlignmentAAList[$endIndex] eq '-' || $variantAlignmentAAList[$endIndex] eq '-'));
			$clusterAA = "$clusterAA*" if($terminated);
			$variantAA = "$variantAA*" if($terminated);
			while($clusterAA ne '' && $variantAA ne '' && isMatched(substr($clusterAA, 0, 1), substr($variantAA, 0, 1))) {
				substr($clusterAA, 0, 1, '');
				substr($variantAA, 0, 1, '');
				$clusterPosition = $clusterPosition + 1;
				$variantPosition = $variantPosition + 1;
			}
			while($clusterAA ne '' && $variantAA ne '' && isMatched(substr($clusterAA, -1, 1), substr($variantAA, -1, 1))) {
				substr($clusterAA, -1, 1, '');
				substr($variantAA, -1, 1, '');
			}
			next if($clusterAA eq $variantAA);
			next if($clusterPosition == 1 && $variantPosition == 1 && length($clusterAA) == 1 && length($variantAA) == 1);
			push(@genotypeList, [join('|', $cluster, $clusterPosition, $clusterAA, $variantAA), $variantPosition]);
		}
		if(@genotypeList) {
			print $writer join("\t", join('|', $tokenHash{'qname'}, $_->[1]), $_->[0]), "\n" foreach(@genotypeList);
		} else {
			print $writer join("\t", $tokenHash{'qname'}, $cluster), "\n";
		}
	}
	close($reader);
	close($writer);

	sub isMatched {
		my ($clusterAA, $variantAA) = @_;
		return 0 if($clusterAA eq '*' || $variantAA eq '*');
		return 1 if($clusterAA eq $variantAA);
		if($allVariants eq '') {
			return 1 if(defined($_ = $aaTypeHash{$variantAA}) && $clusterAA eq $_);
			return 1 if($clusterAA eq '.' && $variantAA ne '-');
			return 1 if($clusterAA eq '_');
		}
		return 0;
	}
}
if($samFile eq '') {
	system("rm $_") foreach("$temporaryPrefix.sam", "$temporaryPrefix.filtered.sam", "$temporaryPrefix.bestfit.sam");
} else {
	system("rm $_") foreach("$temporaryPrefix.sam", "$temporaryPrefix.filtered.sam");
	system("mv $temporaryPrefix.bestfit.sam $samFile");
}
if($alignmentFile ne '') {
	open(my $writer, "| sort -u > $alignmentFile");
	print $writer join("\t", @$_), "\n" foreach(@alignmentTokenListList);
	close($writer);
}

sub getPositionList {
	my ($position, $cigar) = @_;
	my @positionList = ();
	my $index = 0;
	while($cigar =~ s/^([0-9]+)([MIDNSHP=X])//) {
		my ($length, $operation) = ($1, $2);
		if($operation eq 'M') {
			@positionList[$index .. $index + $length - 1] = $position .. $position + $length - 1;
			$index += $length;
			$position += $length;
		} elsif($operation eq 'I') {
			$index += $length;
		} elsif($operation eq 'D') {
			$position += $length;
		} elsif($operation eq 'N') {
			$position += $length;
		} elsif($operation eq 'S') {
			$index += $length;
		}
	}
	return @positionList;
}
