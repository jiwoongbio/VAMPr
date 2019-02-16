# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Cwd 'abs_path';
use Getopt::Long qw(:config no_ignore_case);
use List::Util qw(all);
use XML::LibXML;
use IPC::Open2;
use Bio::DB::Fasta;

(my $vampPath = abs_path($0)) =~ s/\/[^\/]*$//;
my $dataPath = "$vampPath/VAMP_data";

GetOptions('h' => \(my $help = ''),
	'r' => \(my $redownload = ''),
);
if($help) {
	die <<EOF;

Usage:   perl VAMP_database.feature.pl [options]

Options: -h       display this help message
         -r       redownload data

EOF
}

my %proteinClusterHash = ();
{
	open(my $reader, "$dataPath/protein.fasta");
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ s/^>//) {
			next if(scalar(my ($protein, $cluster, $clusterStart, $clusterEnd, $clusterCigar) = split(/\|/, $line)) < 5);
			$proteinClusterHash{$protein} = [$cluster, $clusterStart, $clusterEnd, $clusterCigar];
		}
	}
	close($reader);
}

{
	my $db = Bio::DB::Fasta->new("$dataPath/cluster.fasta");
	sub getClusterSequence {
		my ($cluster, $start, $end) = @_;
		return $db->seq($cluster, $start, $end);
	}
}
my $pid = open2(my $reader, my $writer, "sort -t'\t' -k1,1 -k2,2n -k3,3n -k4,4n -k5,5 -k6,6 -k7,7 -k8,8 -k9,9 | uniq");
foreach my $section ('uniprot_sprot', 'uniprot_trembl') {
	my $URL = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/$section.xml.gz";
	my $file = "$dataPath/$section.xml.gz";
	system("wget --no-verbose -O $file $URL") if(not -r $file or $redownload);

	my $xmlString = '';
	open(my $reader, "gzip -dc $file |");
	while(my $line = <$reader>) {
		if($line =~ /^<entry/) {
			$xmlString = $line;
		} else {
			$xmlString .= $line;
		}
		if($line =~ /^<\/entry/) {
			my @clusterList = ();
			my @proteinList = ();
			push(@proteinList, $1) while($xmlString =~ /<accession>(.*)<\/accession>/g);
			if(my @clusterList = grep {defined} @proteinClusterHash{@proteinList}) {
				my @featureList = getFeatureList($xmlString);
				foreach(@clusterList) {
					my ($cluster, $clusterStart, $clusterEnd, $clusterCigar) = @$_;
					my @clusterPositionList = getPositionList($clusterStart, $clusterCigar);
					foreach(@featureList) {
						my ($start, $end, $original, $variation, $type, $description) = @$_;
						($start, $end) = @clusterPositionList[$start - 1, $end - 1];
						my $clusterSequence = '';
						$clusterSequence = getClusterSequence($cluster, $start, $end) if($original ne '' || $type eq 'sequence variant' || $type =~ /site$/ || $type =~ /residue$/);
						print $writer join("\t", split(/\./, $cluster), $start, $end, $clusterSequence, $original, $_, $type, $description), "\n" foreach($variation ne '' ? split(/,/, $variation, -1) : $variation);
					}
				}
			}
		}
	}
	close($reader);
}
close($writer);
{
	open(my $writer, "> $dataPath/cluster.feature.txt");
	while(my $line = <$reader>) {
		chomp($line);
		$line =~ s/\t/./;
		print $writer "$line\n";
	}
	close($writer);
}
close($reader);
waitpid($pid, 0);

sub getFeatureList {
	my ($xmlString) = @_;
	my @featureList = ();
	my $dom = XML::LibXML->load_xml(string => $xmlString);
	my $root = $dom->documentElement();
	foreach my $featureNode (getChildNodeList($root, 'feature')) {
		my ($type, $description) = map {defined($_) ? $_ : ''} map {$featureNode->getAttribute($_)} ('type', 'description');
		my @positionList = map {$_->getAttribute('position')} getChildNodeList($featureNode, 'location', 'position');
		my @beginList = map {$_->getAttribute('position')} getChildNodeList($featureNode, 'location', 'begin');
		my @endList = map {$_->getAttribute('position')} getChildNodeList($featureNode, 'location', 'end');
		my $original = join(',', map {$_->textContent} getChildNodeList($featureNode, 'original'));
		my $variation = join(',', map {$_->textContent} getChildNodeList($featureNode, 'variation'));
		push(@featureList, [@positionList, @positionList, @beginList, @endList, $original, $variation, $type, $description]);
	}
	return grep {all {defined} @$_} grep {scalar(@$_) == 6} @featureList;
}

sub getChildNodeList {
	my ($node, @childNodeTagNameList) = @_;
	my @childNodeList = ();
	foreach my $childNode ($node->getChildrenByTagName(shift @childNodeTagNameList)) {
		if(@childNodeTagNameList) {
			push(@childNodeList, getChildNodeList($childNode, @childNodeTagNameList));
		} else {
			push(@childNodeList, $childNode);
		}
	}
	return @childNodeList;
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
