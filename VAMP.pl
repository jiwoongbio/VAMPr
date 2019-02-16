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

my @codonList = ();
GetOptions('h' => \(my $help = ''),
	't=s' => \(my $temporaryDirectory = defined($ENV{'TMPDIR'}) ? $ENV{'TMPDIR'} : '/tmp'),
	'C=s' => \@codonList,
	'S=s' => \(my $startCodons = 'GTG,ATG,CTG,TTG,ATA,ATC,ATT'),
	'T=s' => \(my $terminationCodons = 'TAG,TAA,TGA'),
	'L=i' => \(my $minimumTranslationLength = 10),
	'p=i' => \(my $threads = 1),
	'e=f' => \(my $evalue = 10),
	'c=f' => \(my $minimumCoverage = 0.8),
	's=s' => \(my $samFile = ''),
	'a=s' => \(my $alignmentFile = ''),
	'A' => \(my $allVariants = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl VAMP.pl [options] genome.fasta [variant.vcf [sample [...]]] > VAMP.txt

Options: -h       display this help message
         -t DIR   directory for temporary files [\$TMPDIR or /tmp]
         -C STR   codon and translation e.g. ATG=M [NCBI genetic code 11 (Bacterial, Archaeal and Plant Plastid)]
         -S STR   comma-separated start codons [$startCodons]
         -T STR   comma-separated termination codons [$terminationCodons]
         -L INT   minimum translation length [$minimumTranslationLength]
         -p INT   number of threads [$threads]
         -e FLOAT maximum e-value to report alignments [$evalue]
         -c FLOAT minimum coverage [$minimumCoverage]
         -s FILE  output SAM file
         -a FILE  output alignment file
         -A       all variants

EOF
}
my ($inputFastaFile, $inputVcfFile, @inputSampleList) = @ARGV;
my @chromosomeList = ();
my %chromosomeSequenceHash = ();
{
	my $chromosome = '';
	open(my $reader, $inputFastaFile);
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^>(\S*)/) {
			push(@chromosomeList, $chromosome = $1);
		} else {
			$chromosomeSequenceHash{$chromosome} .= uc($line);
		}
	}
	close($reader);
}
chomp(my $hostname = `hostname`);
my $temporaryPrefix = "$temporaryDirectory/VAMP.$hostname.$$";
{
	my %codonHash = (
		'TTT' => 'F', 'CTT' => 'L', 'ATT' => 'I', 'GTT' => 'V',
		'TTC' => 'F', 'CTC' => 'L', 'ATC' => 'I', 'GTC' => 'V',
		'TTA' => 'L', 'CTA' => 'L', 'ATA' => 'I', 'GTA' => 'V',
		'TTG' => 'L', 'CTG' => 'L', 'ATG' => 'M', 'GTG' => 'V',

		'TCT' => 'S', 'CCT' => 'P', 'ACT' => 'T', 'GCT' => 'A',
		'TCC' => 'S', 'CCC' => 'P', 'ACC' => 'T', 'GCC' => 'A',
		'TCA' => 'S', 'CCA' => 'P', 'ACA' => 'T', 'GCA' => 'A',
		'TCG' => 'S', 'CCG' => 'P', 'ACG' => 'T', 'GCG' => 'A',

		'TAT' => 'Y', 'CAT' => 'H', 'AAT' => 'N', 'GAT' => 'D',
		'TAC' => 'Y', 'CAC' => 'H', 'AAC' => 'N', 'GAC' => 'D',
		'TAA' => '*', 'CAA' => 'Q', 'AAA' => 'K', 'GAA' => 'E',
		'TAG' => '*', 'CAG' => 'Q', 'AAG' => 'K', 'GAG' => 'E',

		'TGT' => 'C', 'CGT' => 'R', 'AGT' => 'S', 'GGT' => 'G',
		'TGC' => 'C', 'CGC' => 'R', 'AGC' => 'S', 'GGC' => 'G',
		'TGA' => '*', 'CGA' => 'R', 'AGA' => 'R', 'GGA' => 'G',
		'TGG' => 'W', 'CGG' => 'R', 'AGG' => 'R', 'GGG' => 'G',
	);
	$codonHash{$_->[0]} = $_->[1] foreach(map {[split(/=/, $_)]} @codonList);

	sub translate {
		my ($sequence) = @_;
		return join('', map {defined($_) ? $_ : 'X'} map {$codonHash{substr($sequence, $_ * 3, 3)}} 0 .. int(length($sequence) / 3) - 1);
	}
}
{
	my %startCodonHash = map {$_ => 1} split(/,/, $startCodons);
	my %terminationCodonHash = map {$_ => 1} split(/,/, $terminationCodons);
	my $minimumLength = $minimumTranslationLength * 3;
	open(my $writer, "> $temporaryPrefix.fasta");
	foreach my $chromosome (@chromosomeList) {
		my $sequence = $chromosomeSequenceHash{$chromosome};
		my $sequenceLength = length($sequence);
		foreach my $frame (0 .. 2) {
			my @startIndexList = ();
			for(my $index = $frame; $index + 3 <= $sequenceLength; $index += 3) {
				my $codon = substr($sequence, $index, 3);
				if($startCodonHash{$codon}) {
					push(@startIndexList, $index);
				} elsif(@startIndexList && $terminationCodonHash{$codon}) {
					writeTranslationSequence($chromosome, $sequence, $sequenceLength, '+', $index + 3, @startIndexList);
					@startIndexList = ();
				}
			}
		}
		$sequence = reverseComplementary($sequence);
		foreach my $frame (0 .. 2) {
			my @startIndexList = ();
			for(my $index = $frame; $index + 3 <= $sequenceLength; $index += 3) {
				my $codon = substr($sequence, $index, 3);
				if($startCodonHash{$codon}) {
					push(@startIndexList, $index);
				} elsif(@startIndexList && $terminationCodonHash{$codon}) {
					writeTranslationSequence($chromosome, $sequence, $sequenceLength, '-', $index + 3, @startIndexList);
					@startIndexList = ();
				}
			}
		}
	}
	close($writer);

	sub writeTranslationSequence {
		my ($chromosome, $sequence, $sequenceLength, $strand, $endIndex, @startIndexList) = @_;
		if((my $length = $endIndex - $startIndexList[0]) >= $minimumLength) {
			my ($start, $end) = ('', '');
			($start, $end) = ($startIndexList[0] + 1, $endIndex) if($strand eq '+');
			($start, $end) = (($sequenceLength - $endIndex) + 1, ($sequenceLength - $startIndexList[0])) if($strand eq '-');
			my @startList = map {($_ - $startIndexList[0]) / 3 + 1} @startIndexList;
			print $writer '>', join('|', $chromosome, $start, $end, $strand, @startList), "\n";
			(my $translationSequence = translate(substr($sequence, $startIndexList[0], $length))) =~ s/\*$//;
			print $writer "$translationSequence\n";
		}
	}

	sub reverseComplementary {
		my ($sequence) = @_;
		($sequence = reverse($sequence)) =~ tr/ACGT/TGCA/;
		return $sequence;
	}

	sub terminate {
		my ($sequence) = @_;
		my $sequenceLength = length($sequence);
		for(my $index = 0; $index + 3 <= $sequenceLength; $index += 3) {
			my $codon = substr($sequence, $index, 3);
			return substr($sequence, 0, $index + 3) if($terminationCodonHash{$codon});
		}
		return $sequence;
	}
}
{
#	system("diamond blastp --threads $threads --db $dataPath/protein.dmnd --query $temporaryPrefix.fasta --out $temporaryPrefix.sam --outfmt 101 --max-target-seqs 0 --evalue $evalue --unal 0 --tmpdir $temporaryDirectory --quiet");
	system("diamond blastp --threads $threads --db $dataPath/protein.dmnd --query $temporaryPrefix.fasta --out $temporaryPrefix.sam --outfmt 101 --top 0 --evalue $evalue --unal 0 --tmpdir $temporaryDirectory --quiet");
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
			{
				my ($chromosome, $start, $end, $strand, @startList) = split(/\|/, $tokenHash{'qname'});
				if($tokenHash{'ZS:i'} > 1) {
					my $startIndex = max(grep {$_ <= $tokenHash{'ZS:i'}} @startList) - 1;
					$start += $startIndex * 3 if($strand eq '+');
					$end   -= $startIndex * 3 if($strand eq '-');
					$tokenHash{'ZS:i'} -= $startIndex;
				}
				$tokenHash{'qname'} = join('|', $chromosome, $start, $end, $strand);
			}
			@tagTypeValueList = map {join(':', $_, $tokenHash{$_})} map {join(':', @$_[0, 1])} map {[split(/:/, $_, 3)]} @tagTypeValueList;
			print $writer join("\t", @tokenHash{@samMandatoryFieldList}, @tagTypeValueList), "\n";
		}
	}
	close($reader);
	close($writer);
}
{
	my %variantHash = ();
	if(defined($inputVcfFile)) {
		my $headerLine = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
		my @sampleList = ();
		open(my $reader, $inputVcfFile);
		while(my $line = <$reader>) {
			chomp($line);
			next if($line =~ /^##/);
			next if($line =~ s/^#$headerLine\t// && (@sampleList = split(/\t/, $line)));
			my %tokenHash = ();
			my %sampleTokenHash = ();
			(@tokenHash{split(/\t/, $headerLine)}, @sampleTokenHash{@sampleList}) = split(/\t/, $line);
			if(defined($_ = $tokenHash{'FORMAT'}) && (my @formatList = split(/:/, $_))) {
				foreach my $sample (@sampleList) {
					my %tokenHash = ();
					@tokenHash{@formatList} = split(/:/, $sampleTokenHash{$sample});
					$sampleTokenHash{$sample} = \%tokenHash;
				}
			}
			next unless($tokenHash{'FILTER'} eq 'PASS' || $tokenHash{'FILTER'} eq '.');
			my %genotypeHash = ();
			$genotypeHash{$_} = 1 foreach(map {$_->{'GT'}} @sampleTokenHash{@inputSampleList ? @inputSampleList : @sampleList});
			my @baseList = ($tokenHash{'REF'}, split(/,/, $tokenHash{'ALT'}));
			foreach my $genotype (sort {$a <=> $b} grep {$_ > 0} keys %genotypeHash) {
				my ($chromosome, $position, $refBase, $altBase) = extendIndel(@tokenHash{'CHROM', 'POS'}, @baseList[0, $genotype]);
				my ($start, $end) = ($position, $position + length($refBase) - 1);
				my @indexList = map {int($_ / 1000)} ($start, $end);
				@indexList = @indexList[0, grep {$indexList[$_ - 1] != $indexList[$_]} 1 .. $#indexList] if(scalar(@indexList) > 1);
				my $variant = join("\t", $start, $end, $refBase, $altBase);
				push(@{$variantHash{$chromosome}->{$_}}, $variant) foreach(@indexList);
			}
		}
		close($reader);
	}

	sub getVariantList {
		my ($chromosome, $start, $end, $strand) = @_;
		my @indexList = map {int($_ / 1000)} ($start, $end);
		@indexList = @indexList[0, grep {$indexList[$_ - 1] != $indexList[$_]} 1 .. $#indexList] if(scalar(@indexList) > 1);
		my @variantList = sort map {@$_} grep {defined} map {$variantHash{$chromosome}->{$_}} @indexList;
		@variantList = @variantList[0, grep {$variantList[$_ - 1] ne $variantList[$_]} 1 .. $#variantList] if(scalar(@variantList) > 1);
		my $chromosomeSequence = $chromosomeSequenceHash{$chromosome};
		my @translationVariantList = ();
		foreach(grep {$start <= $_->[0] && $_->[1] <= $end} map {[split(/\t/, $_, -1)]} @variantList) {
			my ($variantStart, $variantEnd, $refBase, $altBase) = @$_;
			if($strand eq '+') {
				my ($translationStart, $translationEnd) = map {int(($_ - $start) / 3) + 1} ($variantStart, (length($refBase) - length($altBase)) % 3 == 0 ? $variantEnd : $end);
				my ($position, $length) = ($start + ($translationStart - 1) * 3, ($translationEnd - $translationStart + 1) * 3);
				my $sequence = (length($refBase) - length($altBase)) % 3 == 0 ? substr($chromosomeSequence, $position - 1, $length) : substr($chromosomeSequence, $position - 1);
				substr(my $variantSequence = $sequence, $variantStart - $position, length($refBase), $altBase);
				my ($translationSequence, $translationVariantSequence) = map {translate(terminate($_))} ($sequence, $variantSequence);
				push(@translationVariantList, [$translationStart, $translationEnd, $translationSequence, $translationVariantSequence]) if($translationSequence ne $translationVariantSequence);
			}
			if($strand eq '-') {
				my ($translationStart, $translationEnd) = map {int(($end - $_) / 3) + 1} ($variantEnd, (length($refBase) - length($altBase)) % 3 == 0 ? $variantStart : $start);
				my ($position, $length) = ($end - ($translationStart - 1) * 3, ($translationEnd - $translationStart + 1) * 3);
				my $sequence = (length($refBase) - length($altBase)) % 3 == 0 ? substr($chromosomeSequence, $position - $length, $length) : substr($chromosomeSequence, 0, $position);
				($sequence, $refBase, $altBase) = map {reverseComplementary($_)} ($sequence, $refBase, $altBase);
				substr(my $variantSequence = $sequence, $position - $variantEnd, length($refBase), $altBase);
				my ($translationSequence, $translationVariantSequence) = map {translate(terminate($_))} ($sequence, $variantSequence);
				push(@translationVariantList, [$translationStart, $translationEnd, $translationSequence, $translationVariantSequence]) if($translationSequence ne $translationVariantSequence);
			}
		}
		return sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @translationVariantList;
	}
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
	{
		sub getVariantSequence {
			my ($chromosome, $start, $end, $strand) = @_;
			my $sequence = substr($chromosomeSequenceHash{$chromosome}, $start - 1, $end - ($start - 1));
			($sequence = reverse($sequence)) =~ tr/ACGT/TGCA/ if($strand eq '-');
			($sequence = translate($sequence)) =~ s/\*$//;
			return $sequence;
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
		my ($chromosome, $start, $end, $strand) = split(/\|/, $tokenHash{'qname'});
		next if(scalar(my ($protein, $cluster, $clusterStart, $clusterEnd, $clusterCigar) = split(/\|/, $tokenHash{'rname'})) < 5);
		my @proteinIndexList = map {defined($_) ? $_ - 1 : $_} getPositionList(@tokenHash{'pos', 'cigar'});
		my @clusterIndexList = map {defined($_) ? $_ - 1 : $_} getPositionList(1, $clusterCigar);
		my @clusterVariantIndexList = grep {defined($_->[0]) && defined($_->[0] = $clusterIndexList[$_->[0]])} map {[$proteinIndexList[$_], $_ + ($tokenHash{'ZS:i'} - 1)]} 0 .. $#proteinIndexList;
		my $clusterSequence = getClusterSequence($cluster, $clusterStart, $clusterEnd);
		my $variantSequence = getVariantSequence($chromosome, $start, $end, $strand);
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
		my %alignmentIndexHash = ();
		foreach(getVariantList($chromosome, $start, $end, $strand)) {
			my $terminated = ($_->[2] =~ /\*$/ || $_->[3] =~ /\*$/);
			my ($startIndex, $endIndex) = ($variantPositionIndexHash{$_->[0]}, $terminated ? $alignmentLength - 1 : $variantPositionIndexHash{$_->[1]});
			my $clusterAA = join('', @clusterAlignmentAAList[$startIndex .. $endIndex]);
			my $variantAA = $_->[3];
			{
				while($startIndex - 1 >= 0 and not isMatched($clusterAlignmentAAList[$startIndex - 1], $variantAlignmentAAList[$startIndex - 1])) {
					$clusterAA = $clusterAlignmentAAList[$startIndex - 1] . $clusterAA;
					$variantAA = $variantAlignmentAAList[$startIndex - 1] . $variantAA;
					$startIndex = $startIndex - 1;
				}
				while($clusterAA ne '' && $variantAA ne '' && isMatched(substr($clusterAA, 0, 1), substr($variantAA, 0, 1))) {
					substr($clusterAA, 0, 1, '');
					substr($variantAA, 0, 1, '');
					$startIndex = $startIndex + 1;
				}
			}
			unless($terminated) {
				while($endIndex + 1 <= $alignmentLength - 1 and not isMatched($clusterAlignmentAAList[$endIndex + 1], $variantAlignmentAAList[$endIndex + 1])) {
					$clusterAA = $clusterAA . $clusterAlignmentAAList[$endIndex + 1];
					$variantAA = $variantAA . $variantAlignmentAAList[$endIndex + 1];
					$endIndex = $endIndex + 1;
				}
				while($clusterAA ne '' && $variantAA ne '' && isMatched(substr($clusterAA, -1, 1), substr($variantAA, -1, 1))) {
					substr($clusterAA, -1, 1, '');
					substr($variantAA, -1, 1, '');
					$endIndex = $endIndex - 1;
				}
				$terminated = ($endIndex == $alignmentLength - 1 && ($clusterAlignmentAAList[$endIndex] eq '-' || $variantAlignmentAAList[$endIndex] eq '-'));
				$variantAA = "$variantAA*" if($terminated);
			}
			my $clusterPosition = $clusterPositionList[$startIndex];
			my $variantPosition = $variantPositionList[$startIndex];
			if($startIndex == $alignmentLength) {
				$clusterPosition = $clusterPositionList[-1] + 1;
				$variantPosition = $variantPositionList[-1] + 1;
			}
			$clusterAA =~ s/-//g;
			$variantAA =~ s/-//g;
			$clusterAA = "$clusterAA*" if($terminated);
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
			$alignmentIndexHash{$_} = 1 foreach($startIndex .. $endIndex);
		}
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
			next if(grep {$alignmentIndexHash{$_}} $startIndex .. $endIndex);
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
			$alignmentIndexHash{$_} = 1 foreach($startIndex .. $endIndex);
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
	system("rm $_") foreach("$temporaryPrefix.fasta", "$temporaryPrefix.sam", "$temporaryPrefix.filtered.sam", "$temporaryPrefix.bestfit.sam");
} else {
	system("rm $_") foreach("$temporaryPrefix.fasta", "$temporaryPrefix.sam", "$temporaryPrefix.filtered.sam");
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

sub extendIndel {
	my ($chromosome, $position, $refBase, $altBase) = @_;
	if($refBase ne $altBase) {
		my $chromosomeSequence = $chromosomeSequenceHash{$chromosome};
		my $isDeletion = $refBase =~ /^$altBase/ || $refBase =~ /$altBase$/;
		while($refBase =~ /^$altBase/ || $altBase =~ /^$refBase/) {
			my $extBase = substr($chromosomeSequence, ($position + length($refBase)) - 1, 1);
			($refBase, $altBase) = map {"$_$extBase"} ($refBase, $altBase);
		}
		while($refBase =~ /$altBase$/ || $altBase =~ /$refBase$/) {
			my $extBase = substr($chromosomeSequence, ($position = $position - 1) - 1, 1);
			($refBase, $altBase) = map {"$extBase$_"} ($refBase, $altBase);
		}
		if($isDeletion) {
			$position += 1;
			substr($_, -1, 1, '') foreach($refBase, $altBase);
			substr($_,  0, 1, '') foreach($refBase, $altBase);
		}
	}
	return ($chromosome, $position, $refBase, $altBase);
}
