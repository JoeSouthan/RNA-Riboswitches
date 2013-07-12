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
use JSON;

use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw (SplitGenBank ProcessSeqs);

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
			print OUT $gb . "\n//\n";
			close (OUT) or die ("Can't close output for $acc\n");
		} else {
			print "Can't get accession!\n";
		}
	}

	print "\nDone\n";
}

#$family, %result
# sub ProcessSequences {
# 	my $family = $_[0];
# 	my %rfam_result = %{$_[1]};

# 	#Check for existing raw files/existing processed files
# 	my $existing_proc = "../output/processed/$family.txt";
# 	my $use_existing_proc = 0;
# 	if (-e $existing_proc){
# 		print "A processed file \"$existing_proc\" exists. Use this?\t";
# 		my $overwrite_choice = <STDIN>;
# 		chomp $overwrite_choice;
# 		if ($overwrite_choice eq "y"){
# 			unlink ($existing_proc);
# 			_processSeqs($family, \%rfam_result);
# 		} else{
# 			$use_existing_proc = 1;
# 		}
# 	} else {
# 		_processSeqs($family, \%rfam_result);
# 	}
# }

#Fam, @gis
sub ProcessSeqs {
	my $family = $_[0];
	my %rfam_result = %{$_[1]};
	my %stats = (
		"RS_Neg" => 0,
		"RS_Pos" => 0,
		"CDS_NoName" => 0
	);
	for my $rfam_fam (keys %rfam_result) {
		#Load the genbank file
		my $fam_mod = $1 if $rfam_fam =~/(.*)\./;
		
		print "\nProcessing $fam_mod\n";
		#BioPerl Object
		eval { 
			my $gb_io = Bio::SeqIO->new(-format => 'genbank', -file => "../output/raw/$family/$fam_mod.txt" ); 
			my $seq_ob = $gb_io->next_seq;
			my $full_Seq = $seq_ob->seq;

			#Counts for loop
			my $count = 0;
			my $cds_count = 0;
			my %gb_cds;
			for my $feature ($seq_ob->get_SeqFeatures) {
				if ($feature->primary_tag eq "CDS") {
					#Get the gene name 
					my $gene_name;
					if ($feature->has_tag("gene")){
						for my $val ($feature->get_tag_values("gene")){
							$gene_name = $val;
						}
					} else {
						$gene_name = $count;
						$count++;
						$stats{"CDS_NoName"}++;
					}
					#Package the features
					my @feat_array = ($feature->start, $feature->end);
					$gb_cds{$gene_name} = \@feat_array;
					$cds_count++;
				}
			}
			$stats{"CDS_Count"} = $cds_count;
			#Just in case there is no CDS in the file
				for my $rs ($rfam_result{$rfam_fam}){
					my $fam_count = 0;
					foreach my $rbs (@{$rs}){
						my $rfam_start = @{$rbs}[0];
						my $rfam_end = @{$rbs}[1];
						my $local_seq = $full_Seq;
						if (0 == $cds_count) {
							#There are no CDS
							print " No features in the GenBank file. Predicting instead...\n";
							my $local_seq = $full_Seq;
							if ("-" eq @{$rbs}[2]){
								$local_seq = _reverseSeq($local_seq);
								$rfam_start = @{$rbs}[1];
								$rfam_end = @{$rbs}[0];
								$stats{"RS_Neg"}++
							} else {
								$stats{"RS_Pos"}++;
							}
							my $proc_seq = _CutSeqPredict($rfam_end, $local_seq, $rfam_start);
							_OutputProc ($family, $proc_seq, $fam_count, $rfam_fam, "");

							my @stat_array = $rfam_fam.$fam_count;
							$stats{"CDS_Fam_NoFeat"} = \@stat_array;
							

						} else {

							#There are CDS
							my $current_lowest = 999999;
							my $chosen_locus;

							#If negative, reverse the numbers and sequence
							if ("-" eq @{$rbs}[2]){
								$local_seq = _reverseSeq($local_seq);
								$rfam_start = @{$rbs}[1];
								$rfam_end = @{$rbs}[0];
								$stats{"RS_Neg"}++;
							} else {
								$stats{"RS_Pos"}++;
							}

							#Go through the GenBank File
							for my $locus (keys (%gb_cds)){
								my $difference = $gb_cds{$locus}[0]-$rfam_start;
								if ($current_lowest > $difference and $difference >= 0 ) {
									$current_lowest = $difference;
									$chosen_locus = $locus;
								}
							}

							# In the case where there are no substantial differences 
							# between the genes, send it for prediction instead
							if (!defined($chosen_locus)) {
								print " No significant differences, predicting...\n";
								my $proc_seq = _CutSeqPredict($rfam_end, $local_seq, $rfam_start);
								_OutputProc ($family, $proc_seq, $fam_count, $rfam_fam, "");
							} else {
								my $length_req = (abs($rfam_end-$rfam_start))+$current_lowest;
								print " Chosen locus $chosen_locus at difference $current_lowest length $length_req rfamlength".abs($rfam_end-$rfam_start)."\n";
								my $proc_seq = _CutSeq($rfam_start, $length_req, $local_seq);
								_OutputProc ($family, $proc_seq, $fam_count, $rfam_fam, $chosen_locus);
							}
						}
						$fam_count++;
					}
				}
		};
		#print "*Can't open sequence for $rfam_fam $@\n" if $@;
		print "*Can't open sequence for $rfam_fam\n" if $@;
		#Error handling? 
	}
	_OutputStats(\%stats, $family);
}
#Reverses a sequence 
sub _reverseSeq {
	my $seq = $_[0];
	$seq =~ tr/acgtrymkbdhvACGTRYMKBDHV/tgcayrkmvhdbTGCAYRKMVHDB/;
	return $seq;
}
sub _CutSeq {
	my ($rfam_start, $length_req, $seq) = @_;
	my $processed_seq = substr ($seq, $rfam_start, $length_req+3 );
	return $processed_seq;
}
sub _CutSeqPredict {
	my $rfam_end = $_[0];
	my $seq = $_[1];
	my $rfam_start = $_[2];

	my $meth_site = index($seq, "ATG", $rfam_end);
	my $seq_length = (abs($rfam_end-$rfam_start))+($rfam_end-$meth_site);
	my $processed_seq = substr ($seq, $rfam_start, $seq_length);
	return $processed_seq;
}
sub _OutputProc {
	my ($family, $processed_seq, $fam_count, $rfam_fam, $locus)  = @_;
	make_path("../output/processed/$family");

	my $filename = $rfam_fam."_".$fam_count;
	open (OUTPUT, ">", "../output/processed/$family/$filename.txt") or die "Can't create output file\n";
	print OUTPUT ">$filename | ".length($processed_seq)."bp | $locus\n";
	print OUTPUT $processed_seq."\n";
	close OUTPUT or die "Can't close output\n";
}
sub _OutputStats {
	my %stats = %{$_[0]};
	my $family = $_[1];
	my $json = JSON->new;
	my $encoded = $json->pretty->encode(\%stats);

	make_path("../output/stats");

	open (STATS, ">", "../output/stats/$family.txt") or die "Can't create Stats\n";
	print STATS $encoded;
	close STATS or die "Can't close stat file\n";
}
1;
