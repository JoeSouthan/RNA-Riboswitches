#! /usr/bin/perl -w
#
# getseqs.pl 
# Created by: Joseph Southan
# Date: 23/5/13
# Description: - Retrieves DNA/RNA sequences from RFAM
# Usage: getseqs.pl [seed/full] [text file]
#				Seed/Full = To get the full aligment file or just the seeds (Default: seed)
#				Text file = File containing RFAM families with each one on a newline
#
use strict;
use LWP::Simple;
use Data::Dumper;
use Bio::SeqIO;
use Bio::DB::EUtilities;
use LWP::UserAgent;
use File::Path qw(make_path remove_tree);


#
#
#	Script Logic
#
#
my ($mode, $limit, $search_file, $rate_limit, @rfam_families);

#Turn off buffer for status messages
$|=1;
#Check for params
die ("\n========================\nError: Please ensure to have the correct arguments\n========================\n") unless (@ARGV);
#Default to seed
unless ($ARGV[0] eq "seed" or $ARGV[0] eq "full") {
	$mode = "seed";
	print "\nNote: Wrong mode for RFAM, defaulting to \"seed\"\n\n";
} else {
	$mode = $ARGV[0];
}
unless (defined($ARGV[1])) {
	die ("No search file found\n");
} else {
	#Open up the search file
	$search_file = $ARGV[1];
	open (SEARCHES, "<", $search_file) or die ("Can't open search file, ensure it exists.\n");
	while (my $line = <SEARCHES>){
		chomp ($line);
		push (@rfam_families, $line);
	}
	close (SEARCHES) or die ("Can't close search file.\n");
}
#Warn the user
if ($mode eq "full"){
	print "\n========================\nAre you sure you want to get the full alignments? \nThis is very data and processor hungry and some sequences may not contain the full alignments. y/n \n";
	my $choice = <STDIN>;
	chomp ($choice);
	if ($choice eq "n") {
		$mode = "seed";
	} 
}


#
#
#	Main loop
#
#
my $rfam_families_count = scalar(@rfam_families);
foreach my $family (@rfam_families) {
	system("cls");
	print "\n========================\nFamily: $family $rfam_families_count remain(s).\n========================\n";
	#
	#
	#	RFAM
	#
	#
	#Read in the RFAM list
	#http://rfam.sanger.ac.uk/family/RF00174/alignment/seed/stockholm?alnType=seed&nseLabels=0&cs=0
	print "\n\nRetrieving $mode alignments from RFAM\n========================\n";
	my %rfam_result;
	my $count = 0;
	my $exist_count = 0;
	my $rfam_url = get ("http://rfam.sanger.ac.uk/family/$family/alignment/$mode/stockholm?alnType=$mode&nseLabels=0&cs=0");
	if ($rfam_url) {
		##=GS C.diphtheriae.1         AC    BX248356.1/172244-172444
		while ($rfam_url =~ /([A-Z0-9.]+)\/(\d+)-(\d+)\n/g){
			my @locs = ($2, $3);
			my $polarity;
			if ($3 < $2){
				$polarity = "-";
			} else {
				$polarity = "+";
			}
			#Duplicate key fix
			push ( @locs, $polarity);
			if ($rfam_result{$1}){
				#print "$1 exists already, adding to to existing key";
				my @sub_result = [ @locs ];
				push( @{ $rfam_result{$1}}, @sub_result); 
				$exist_count++;
			} else {
				my @sub_result = [@locs];
				$rfam_result{$1} = [@sub_result];
			}
			$count++;
		}
	} else {
		print "$family Get Failed\n";
		sleep 10;
		$rfam_families_count--;
		next;
	}
	my @gis = keys (%rfam_result);
 	my $no_ids = scalar(@gis);
	print "Found $count records from RFAM for $family family. ($no_ids individual genes)\n";

	

	#
	#
	#	NCBI
	#
	#
	#NCBI FASTA eghttp://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=BX248356.1&rettype=fasta
	#Get the full sequence from the start location
	print "\n\nRetrieving sequences from NCBI for $family \n========================\n\r";
	my $existing_raw = "../output/raw/$family.txt";
	my $existing_raw_split ="../output/raw/$family";
	my $use_existing_raw = 0;
	my $fetch_result = 0;
	if (-e $existing_raw or -d $existing_raw_split) {
		print "Raw files for this family exist. Overwrite?\t";
		my $overwrite_choice = <STDIN>;
		chomp $overwrite_choice;
		if ($overwrite_choice eq "y") {
			unlink($existing_raw);
			if (-d $existing_raw_split) {
				remove_tree($existing_raw_split);
			}
			$fetch_result = FetchSeqs(\@gis, $family);
		} else {
			$use_existing_raw = 1; #
		}
	} else {
		$fetch_result = FetchSeqs(\@gis, $family);
	}


	if (0 == $use_existing_raw) {
		print "\n\nSplitting GenBank files \n========================\n\r";
		SplitGenBank($family, $no_ids);
	}

	#
	#
	# Process the sequences
	#
	#
	#Check for existing raw files/existing processed files
	my $existing_proc = "../output/processed/$family.txt";
	my $use_existing_proc = 0;
	if (-e $existing_proc){
		print "A processed file \"$existing_proc\" exists. Use this?\t";
		my $overwrite_choice = <STDIN>;
		chomp $overwrite_choice;
		if ($overwrite_choice eq "y"){
			unlink ($existing_proc);
			processSeqs($family, \%rfam_result);
		} else{
			$use_existing_proc = 1;
		}
	} else {
		processSeqs($family, \%rfam_result);
	}


	$rfam_families_count--;
	sleep 3;
}
#Reverses a sequence 
sub reverseSeq {
	my $seq = $_[0];
	$seq =~ tr/acgtrymkbdhvACGTRYMKBDHV/tgcayrkmvhdbTGCAYRKMVHDB/;
	return $seq;
}
sub FetchSeqs {
	my @gis = @{$_[0]};
 	my $family = $_[1];
 	my $ua = LWP::UserAgent->new();
 	my $fetch_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gbwithparts&id=";

 	$fetch_url .= join (",", @gis);

 	#Progres bar
 	#$ua->show_progress('true value
 	print "Getting sequences, please wait...\r";
  	my $start_time = time();
 	my $response = $ua->get($fetch_url,':content_file' => "../output/raw/$family.txt");
 	if ($response->is_success){
 		print "Sequences obtained in @{[time()-$start_time]} seconds\n";
 		return 1;
 	} else {
 		print "Error Retrieving sequences: ".$response->status_list."\n";
 		return 0;
 	}
}
#family, no_ids
sub SplitGenBank {
	my ($family, $no_ids) = @_;
	my $raw_file = Bio::SeqIO->new(-format => 'genbank', -file => "../output/raw/$family.txt");
	#Create directory
	make_path("../output/raw/$family");

	#Start splitting the file
	while (my $seq = $raw_file->next_seq){
		my $accession = $seq->accession_number;
		print "Processing $accession\r";
		my $output = Bio::SeqIO->new(-format => 'genbank', -file => ">../output/raw/$family/$accession.txt", -verbose => -1);
		$output->write_seq($seq);
	}
	print "Done\n";
}
#Fam, @gis
sub processSeqs {
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
					##$ncbi_result{$keys} = reverseSeq($ncbi_result{$keys});
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