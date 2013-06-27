#! /usr/bin/perl -w
#
# ProcessSeqs.pm
# Created by: Joseph Southan
# Date: 27/6/13
# Description: - Module for getseqs.pl
# 
package ProcessSeqs;
use strict;
use Data::Dumper;
use Bio::SeqIO;
use File::Path qw(make_path remove_tree);
use Tie::File;

use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw (SplitGenBank ProcessSequences);

#family
sub SplitGenBank {
	my $family = $_[0];

	#Create directory
	make_path("../output/raw/$family");

	tie my @genbank_file, 'Tie::File', "../output/raw/$family.txt", recsep => "//\n";
	foreach my $gb (@genbank_file) {
		my $acc = ($gb =~ /LOCUS\s+([A-Z0-9.]+)\s/g) ? $1 : 0;
		unless (0 eq $acc) {
			print "Outputting: $acc\r";
			open (OUT, ">", "../output/raw/$family/$acc.txt") or die ("Can't open output for $acc\n");
			print OUT $gb . "//\n";
			close (OUT) or die ("Can't close output for $acc\n");
		} else {
			print "Can't get accession!\n";
		}
	}

	print "\nDone\n";
}

#$family, %result
sub ProcessSequences {
	my $family = $_[0];
	my %rfam_result = %{$_[1]};

	#Check for existing raw files/existing processed files
	my $existing_proc = "../output/processed/$family.txt";
	my $use_existing_proc = 0;
	if (-e $existing_proc){
		print "A processed file \"$existing_proc\" exists. Use this?\t";
		my $overwrite_choice = <STDIN>;
		chomp $overwrite_choice;
		if ($overwrite_choice eq "y"){
			unlink ($existing_proc);
			_processSeqs($family, \%rfam_result);
		} else{
			$use_existing_proc = 1;
		}
	} else {
		_processSeqs($family, \%rfam_result);
	}
}

#Fam, @gis
sub _processSeqs {
	my $family = $_[0];
	my %rfam_result = %{$_[1]};
	my $processed_count = 0;
	#Check raw exists
	my $filename = "../output/raw/$family.txt";
	my $filename_out = "../output/processed/$family.txt";
	unless (-e $filename){
		return "$filename does not exist!\n";
	}
	# tie my @raw_file, $filename or die "Can't tie $filename\n";
	#Change delimiter
	local $/ = '>';
	#Read raw file
	open (RAW, "<", $filename) or die "Can't open $filename\n";
	my %seqs;
	while (my $line = <RAW>){
		my @lines = split "\n", $line;
		my $info_line = shift(@lines);
		if ($info_line =~/gi\|(\d+)\|(.*)/){
			print "$1\n";
			my $gi = $1;
			$seqs{$gi} = join ('', @lines);
			#pop (@lines);
		}
		
	}
	close (RAW) or die "Can't close $filename\n";
	open (OUTPUT, ">>", $filename_out);
	print Dumper $seqs{"AAYA01000006.1"};
	die;
	foreach my $fams (keys(%rfam_result)) {
		#print "fams: $fams \n";
		#Get the sequence
		my $sequence = $seqs{$fams};
		print "$sequence $fams";
		die;
		my @fam_arrays = $rfam_result{$fams};
		# $fam_arrays[$i][0][2];
		for (my $i = 0; $i < @fam_arrays; $i++) {
			my @fam_sub_arrays = $fam_arrays[$i];
			for (my $j = 0; $j < @fam_sub_arrays; $j++) {
				my $start = $fam_arrays[$i][$j][0];
				my $end = $fam_arrays[$i][$j][1];
				my $pol = $fam_arrays[$i][$j][2];
				#print "\t$start $end $pol\n";
				if ($pol eq "-") {
					#reverse and transwhatsit
					##$ncbi_result{$keys} = _reverseSeq($ncbi_result{$keys});
					#swap the start and stop
					$start = $fam_arrays[$i][$j][1];
					$end = $fam_arrays[$i][$j][0];
				}
				my $seq_length = $end-$start;
				my $meth_site = index($sequence, "ATG", $end);

				#print "pol: $polarity rs: $rfam_start, re: $rfam_end, SeqL: $seq_length, meth: $meth_site \n";    
				my $processed_seq = substr ($sequence, $start, ($seq_length+($meth_site-$seq_length)+3));
				print OUTPUT $processed_seq;
				#$process_count++;
			}
		}	
	}
	close (OUTPUT) or die "Can't close $filename_out\n";
	# print "Files processed, remove raw sequence file?\t";
	# my $choice = <STDIN>;
	# chomp $choice;
	# if ($choice eq "y") {
	# 	unlink ("../output/raw/$family.txt");
	# 	print "\nFile deleted.\n";
	# } else {
	# 	print "Processed files avaliable in ../output/processed/$family.txt\n";
	# 	last;
	# }
}
#Reverses a sequence 
sub _reverseSeq {
	my $seq = $_[0];
	$seq =~ tr/acgtrymkbdhvACGTRYMKBDHV/tgcayrkmvhdbTGCAYRKMVHDB/;
	return $seq;
}
1;
