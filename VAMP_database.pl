# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Cwd 'abs_path';
use Getopt::Long qw(:config no_ignore_case);
use List::Util qw(max);
use IPC::Open2;
use Bio::DB::Fasta;

my %aaTypeHash = (
	'R' => 'b', 'H' => 'b', 'K' => 'b', # basic (b)
	'D' => 'a', 'E' => 'a', # acidic (a)
	'S' => 'p', 'T' => 'p', 'N' => 'p', 'Q' => 'p', # polar (p)
	'A' => 'h', 'V' => 'h', 'I' => 'h', 'L' => 'h', 'M' => 'h', 'F' => 'h', 'Y' => 'h', 'W' => 'h', # hydrophobic (h)
);

(my $vampPath = abs_path($0)) =~ s/\/[^\/]*$//;
my $dataPath = "$vampPath/VAMP_data";
system("mkdir -p $dataPath");

my $vampURL = 'http://share.zhanxw.com/VAMPr';
my $dataURL = "$vampURL/VAMP_data";

my @additionalOrthologyProteinSequenceFileList = ();
GetOptions('h' => \(my $help = ''),
	'r' => \(my $redownload = ''),
	'c=f' => \(my $sequenceIdentityThreshold = 0.7),
	'd' => \(my $decoy = ''),
	'a=s' => \@additionalOrthologyProteinSequenceFileList,
	'p=i' => \(my $threads = 1),
);
if($help) {
	die <<EOF;

Usage:   perl VAMP_database.pl [options] [KEGG_orthology.txt]

Options: -h       display this help message
         -r       redownload data
         -c FLOAT sequence identity threshold [$sequenceIdentityThreshold]
         -d       decoy
         -a FILE  additional protein sequence file including orthology, protein and sequence columns
         -p INT   number of threads [$threads]

EOF
}
my ($orthologyFile) = @ARGV;
if(defined($orthologyFile)) {
	my %orthologyHash = ();
	{
		open(my $reader, $orthologyFile);
		while(my $line = <$reader>) {
			chomp($line);
			my ($orthology) = split(/\t/, $line);
			$orthologyHash{$orthology} = 1;
		}
		close($reader);
	}

	my $file = $decoy eq '' ? "$dataPath/protein_sequence.txt" : "$dataPath/protein_sequence.decoy.txt";
	writeProteinSequence($file) if(not -r $file or $redownload);

	my %pidFileHash = ();
	open(my $reader, $file);
	open(my $writerProteinAlignment, ">> $dataPath/protein_alignment.txt");
	open(my $writerProtein, ">> $dataPath/protein.fasta");
	open(my $writerCluster, ">> $dataPath/cluster.fasta");
	my ($orthology, %proteinSequenceHash) = ('');
	while(my $line = <$reader>) {
		chomp($line);
		my @tokenList = split(/\t/, $line);
		if($tokenList[0] ne $orthology) {
			writeProteinClusterFasta() if($orthologyHash{$orthology});
			($orthology, %proteinSequenceHash) = ($tokenList[0]);
		}
		$proteinSequenceHash{$tokenList[1]} = $tokenList[2];
	}
	writeProteinClusterFasta() if($orthologyHash{$orthology});
	while((my $pid = wait()) != -1) {
		my $file = $pidFileHash{$pid};
		readProteinSequenceAlignment($file);
		delete $pidFileHash{$pid};
	}
	close($reader);
	close($writerProteinAlignment);
	close($writerProtein);
	close($writerCluster);

	sub writeProteinClusterFasta {
		if($decoy eq '') {
			my $file = "$dataPath/$orthology.protein_alignment.txt";
			if(my $pid = fork()) {
				$pidFileHash{$pid} = $file;
			} else {
				writeProteinSequenceAlignment($file);
				exit(0);
			}
			if(scalar(keys %pidFileHash) == $threads) {
				my $pid = wait();
				my $file = $pidFileHash{$pid};
				readProteinSequenceAlignment($file);
				delete $pidFileHash{$pid};
			}
		} else {
			foreach my $protein (sort keys %proteinSequenceHash) {
				print $writerProtein ">$protein|$orthology\n";
				print $writerProtein "$proteinSequenceHash{$protein}\n";
			}
		}
	}

	sub writeProteinSequenceAlignment {
		my %clusterNumberProteinListHash = ();
		{
			open(my $writer, "> $dataPath/$orthology.fasta");
			print $writer ">$_\n$proteinSequenceHash{$_}\n" foreach(sort keys %proteinSequenceHash);
			close($writer);

			system("cd-hit -i $dataPath/$orthology.fasta -o $dataPath/$orthology.cd-hit -c $sequenceIdentityThreshold -d 0 > /dev/null");

			my $clusterNumber = '';
			open(my $reader, "$dataPath/$orthology.cd-hit.clstr");
			while(my $line = <$reader>) {
				chomp($line);
				if($line =~ /^>Cluster ([0-9]+)$/) {
					$clusterNumber = $1;
				} elsif($line =~ /^[0-9]+\t[0-9]+aa, >(.*)\.\.\. .*$/) {
					push(@{$clusterNumberProteinListHash{$clusterNumber}}, $1);
				} else {
					print STDERR "$line\n";
				}
			}
			close($reader);

			system("rm $dataPath/$orthology.fasta $dataPath/$orthology.cd-hit $dataPath/$orthology.cd-hit.clstr");
		}

		my ($file) = @_;
		open(my $writer, "> $file");
		foreach my $clusterNumber (sort {$a <=> $b} keys %clusterNumberProteinListHash) {
			my $cluster = "$orthology.$clusterNumber";
			my @proteinList = @{$clusterNumberProteinListHash{$clusterNumber}};
			my %proteinAlignmentHash = ();
			if(scalar(@proteinList) > 1) {
#				my $pid = open2(my $reader, my $writer, 'clustalo --infile=- --seqtype=Protein --threads=1');
				my $pid = open2(my $reader, my $writer, 'mafft --auto --quiet --thread 1 --anysymbol /dev/stdin');
				foreach my $protein (@proteinList) {
					print $writer ">$protein\n";
					print $writer "$proteinSequenceHash{$protein}\n";
				}
				close($writer);
				my $protein = '';
				while(my $line = <$reader>) {
					chomp($line);
					if($line =~ s/^>//) {
						($protein = $line) =~ s/ .*$//;
					} else {
						$proteinAlignmentHash{$protein} .= $line;
					}
				}
				close($reader);
				waitpid($pid, 0);
			} else {
				@proteinAlignmentHash{@proteinList} = @proteinSequenceHash{@proteinList};
			}
			print $writer join("\t", $cluster, $_, $proteinSequenceHash{$_}, $proteinAlignmentHash{$_}), "\n" foreach(@proteinList);
		}
		close($writer);
	}

	sub readProteinSequenceAlignment {
		my ($file) = @_;
		open(my $reader, $file);
		my ($cluster, @proteinSequenceAlignmentList) = ('');
		while(my $line = <$reader>) {
			chomp($line);
			my @tokenList = split(/\t/, $line);
			if($tokenList[0] ne $cluster) {
				writeProteinClusterFastaPart($cluster, @proteinSequenceAlignmentList) if($cluster ne '');
				($cluster, @proteinSequenceAlignmentList) = ($tokenList[0]);
			}
			push(@proteinSequenceAlignmentList, [@tokenList[1 .. 3]]);
		}
		writeProteinClusterFastaPart($cluster, @proteinSequenceAlignmentList) if($cluster ne '');
		close($reader);
		system("rm $file");
	}

	sub writeProteinClusterFastaPart {
		my ($cluster, @proteinSequenceAlignmentList) = @_;
		my @alignmentList = ();
		foreach(@proteinSequenceAlignmentList) {
			my ($protein, $sequence, $alignment) = @$_;
			$alignment =~ s/^(-*)//; my $start = length($1) + 1;
			$alignment = (' ' x length($1)) . $alignment;
			$alignment =~ s/(-*)$//; my $end = length($alignment);
			$alignment = $alignment . (' ' x length($1));
			push(@alignmentList, $alignment);
			print $writerProteinAlignment join("\t", $cluster, $protein, $alignment), "\n";
			my $cigar = '';
			while($alignment ne '') {
				if($alignment =~ s/^([A-Z]+)//) {
					$cigar .= length($1) . 'M';
				} elsif($alignment =~ s/^(-+)//) {
					$cigar .= length($1) . 'D';
				} elsif($alignment =~ s/^( +)//) {
				}
			}
			print $writerProtein ">$protein|$cluster|$start|$end|$cigar\n";
			print $writerProtein "$sequence\n";
		}
		my $length = max(map {length($_)} @alignmentList);
		my $clusterSequence = '';
		foreach my $index (0 .. $length - 1) {
			my $aa = '.';
			my @aaList = grep {$_ ne ' '} map {substr($_, $index, 1)} @alignmentList;
			if(scalar(grep {$_ ne $aaList[0]} @aaList) == 0) {
				$aa = $aaList[0];
			} else {
				my @aaTypeList = map {defined($aaTypeHash{$_}) ? $aaTypeHash{$_} : $_} @aaList;
				$aa = $aaTypeList[0] if(scalar(grep {$_ ne $aaTypeList[0]} @aaTypeList) == 0);
			}
			$aa = '_' if(grep {$_ eq '-'} @aaList);
			$clusterSequence .= $aa;
		}
		print $writerCluster ">$cluster\n";
		print $writerCluster "$clusterSequence\n";
	}

	sub writeProteinSequence {
		my %uniprotOrthologyHash = ();
		{
			my $file = $decoy eq '' ? "$dataPath/orthology_uniprot.txt" : "$dataPath/orthology_uniprot.decoy.txt";
			if(not -r $file or $redownload) {
				open(my $writer, "> $file");
				foreach my $orthology (sort keys %orthologyHash) {
					open(my $reader, "wget --no-verbose -O - http://rest.genome.jp/link/uniprot/$orthology |");
					while(my $line = <$reader>) {
						chomp($line);
						my ($orthology, $uniprot) = split(/\t/, $line);
						$orthology =~ s/^ko://;
						$uniprot =~ s/^up://;
						print $writer join("\t", $orthology, $uniprot), "\n";
					}
					close($reader);
				}
				close($writer);
			}
			open(my $reader, $file);
			while(my $line = <$reader>) {
				chomp($line);
				my ($orthology, $uniprot) = split(/\t/, $line);
				$uniprotOrthologyHash{$uniprot}->{$orthology} = 1 if($orthologyHash{$orthology});
			}
			close($reader);
		}
		my %geneOrthologyHash = ();
		{
			my $file = $decoy eq '' ? "$dataPath/orthology_gene.txt" : "$dataPath/orthology_gene.decoy.txt";
			if(not -r $file or $redownload) {
				open(my $writer, "> $file");
				foreach my $orthology (sort keys %orthologyHash) {
					open(my $reader, "wget --no-verbose -O - http://rest.kegg.jp/link/genes/$orthology |");
					while(my $line = <$reader>) {
						chomp($line);
						my ($orthology, $gene) = split(/\t/, $line);
						$orthology =~ s/^ko://;
						print $writer join("\t", $orthology, $gene), "\n";
					}
					close($reader);
				}
				close($writer);
			}
			open(my $reader, $file);
			while(my $line = <$reader>) {
				chomp($line);
				my ($orthology, $gene) = split(/\t/, $line);
				$geneOrthologyHash{$gene}->{$orthology} = 1 if($orthologyHash{$orthology});
			}
			close($reader);
		}
		{
			my $URL = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz";
			my $file = "$dataPath/idmapping.dat.gz";
			system("wget --no-verbose -O $file $URL") if(not -r $file or $redownload);

			my %tokenHash = ();
			open(my $reader, "gzip -dc $file |");
			while(my $line = <$reader>) {
				chomp($line);
				my @tokenList = split(/\t/, $line);
				if($tokenList[1] eq 'KEGG' && defined(my $orthologyHash = $geneOrthologyHash{$tokenList[2]})) {
					my $uniprot = $tokenList[0];
					$uniprotOrthologyHash{$uniprot}->{$_} = 1 foreach(keys %$orthologyHash);
				}
			}
			close($reader);
		}

		my ($file) = @_;
		open(my $writer, "| sort -t'\t' -k1,1 -k2,2 -k3,3 | uniq > $file");
		foreach my $section ('uniprot_sprot', 'uniprot_trembl') {
			my $URL = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/$section.fasta.gz";
			my $file = "$dataPath/$section.fasta.gz";
			system("wget --no-verbose -O $file $URL") if(not -r $file or $redownload);

			my ($uniprot, $sequence) = ('', '');
			open(my $reader, "gzip -dc $file |");
			while(my $line = <$reader>) {
				chomp($line);
				if($line =~ s/^>//) {
					if(defined(my $orthologyHash = $uniprotOrthologyHash{$uniprot})) {
						print $writer join("\t", $_, $uniprot, $sequence), "\n" foreach(keys %$orthologyHash);
					}
					$uniprot = $line;
					$uniprot =~ s/^(sp|tr)[|]//;
					$uniprot =~ s/[| ].*$//;
					$sequence = '';
				} else {
					$sequence .= $line;
				}
			}
			if(defined(my $orthologyHash = $uniprotOrthologyHash{$uniprot})) {
				print $writer join("\t", $_, $uniprot, $sequence), "\n" foreach(keys %$orthologyHash);
			}
			close($reader);
		}
		foreach my $file (@additionalOrthologyProteinSequenceFileList) {
			open(my $reader, $file);
			print $writer $_ while(<$reader>);
			close($reader);
		}
		close($writer);
	}
} else {
	{
		my $URL = "$dataURL/protein.fasta";
		my $file = "$dataPath/protein.fasta";
		system("wget --no-verbose -O $file $URL") if(not -r $file or $redownload);
	}
	{
		my $URL = "$dataURL/cluster.fasta";
		my $file = "$dataPath/cluster.fasta";
		system("wget --no-verbose -O $file $URL") if(not -r $file or $redownload);
	}
}
system("diamond makedb --db $dataPath/protein.dmnd --in $dataPath/protein.fasta");
system("rm -f $dataPath/cluster.fasta.index");
Bio::DB::Fasta->new("$dataPath/cluster.fasta");
