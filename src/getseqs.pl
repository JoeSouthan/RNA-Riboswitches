#! /usr/bin/perl -w
#
# getseqs.pl 
# Created by: Joseph Southan
# Date: 23/5/13
# Description: - Retrieves DNA/RNA sequences from RFAM
# Usage: getseqs.pl [seed/full] [text file] [Rate limit] (Limit)
#				Seed/Full = To get the full aligment file or just the seeds (Default: seed)
#				Text file = File containing RFAM families with each one on a newline
#				Rate limit= Limit the script to fetch x number of sequences at a time (Default: 10 seconds)
#				limit 	  = Limit to only x sequences per family (optional)
#
use strict;
use LWP;
use LWP::Simple;
use Data::Dumper;
use Bio::DB::EUtilities;
use Term::ProgressBar;


#
#
#	Script Logic
#
#
my ($mode, $limit, $search_file, $rate_limit, @rfam_families);

#Turn off buffer for status messages
$|=1;
#Check for params
unless (@ARGV) {
	die ("\n========================\nError: Please ensure to have the correct arguments\n========================\n");
}
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
#Rate limiting to keep us friendly with API providers
unless (defined($ARGV[2])) {
	$rate_limit = 10;
} else {
	$rate_limit = $ARGV[2];
}
#If only x sequences per family are required set a limit
if (defined($ARGV[3])) {
	$limit = $ARGV[3];
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
	unless (defined($limit)) {
		print "\n========================\nFamily: $family $rfam_families_count remain(s).\nSequences at a time: $rate_limit\n========================\n";
	} else {
		print "\n========================\nFamily: $family $rfam_families_count remain(s).\nSequences at a time: $rate_limit Limit: $limit\n========================\n";
	}
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
	print "Found $count Genes from RFAM for $family family.\n";

	#
	#
	# Convert and retrieve NCBI sequences
	# (Adapted from BioPerl documentation)
	#
	my @ids = keys (%rfam_result);
 	print "IDS ".scalar(@ids)."\n";
	my $factory = Bio::DB::EUtilities->new(-eutil   => 'efetch',
                                       -db      => 'nuccore',
                                       -id      => \@ids,
                                       -email   => 'joseph@southanuk.co.uk',
                                       -rettype => 'gi');
 
	my @gis = split(m{\n},$factory->get_Response->content);

	foreach my $genes (@gis) {
		
	}

	#
	#
	#	NCBI
	#
	#
	#NCBI FASTA eghttp://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=BX248356.1&rettype=fasta
	#Get the full sequence from the start location
	print "\n\nRetrieving sequences from NCBI for $family \n========================\n\r";
	my $fetch_batch = 100;
	my $gi_no = scalar(@gis);
	my $g_count = 0;
	for (my $i = 0; $i < $gi_no; $i += $fetch_batch) {
		my @get_ids;
		for (my $j = 0; $j < $fetch_batch; $j++) {
			if (@gis) {
				push (@get_ids, pop(@gis));
				$g_count++;
			} else {
				last;
			}
		}
		my $fetch_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&id=";
		$fetch_url .= join (",", @get_ids);
		
		#die($fetch_url);
		my $get_fetch_url = get($fetch_url);
		if ($get_fetch_url) {
			open (SEQS, ">", "../output/raw/$family.txt") or die "Can't create output file\n";
			print SEQS ">>First 100\n";
			print SEQS $get_fetch_url;
			close SEQS or die "Can't close output file\n";
		} else {
			print "Error getting the sequences\n";
		}
	}
		my (%ncbi_result,@ncbi_result_fail,$duration);
	my $result_count = 0;
	my $fail_count = 0;
	
	#print "\n\nGot result_count sequences in $duration seconds. $fail_count failed attemtps.\n";

	#
	#
	# Process the sequences
	#
	#
	print "\n\nNow Processing sequences for $family...\n========================\n";
	my %processed_results;
	my $process_count = 0;
	foreach my $keys (keys (%ncbi_result)){
		print "$process_count\r";
		my $rfam_start = $rfam_result{$keys}[0];
		my $rfam_end = $rfam_result{$keys}[1];
		my $polarity = $rfam_result{$keys}[2];

		if ($polarity eq "-") {
			#reverse and transwhatsit
			$ncbi_result{$keys} = reverseSeq($ncbi_result{$keys});
			#swap the start and stop
			$rfam_start = $rfam_result{$keys}[1];
			$rfam_end = $rfam_result{$keys}[0];
		}
		my $seq_length = $rfam_end-$rfam_start;
		my $meth_site = index($ncbi_result{$keys}, "ATG", $seq_length);

		#print "pol: $polarity rs: $rfam_start, re: $rfam_end, SeqL: $seq_length, meth: $meth_site \n";
		my $processed_seq = substr ($ncbi_result{$keys}, 0, ($seq_length+($meth_site-$seq_length)+3));
		$processed_results{$keys} = $processed_seq;
		$process_count++;
	}
	print "\r$process_count Sequences processed\n";

	#
	#
	# Output the sequences
	#
	#
	#Output the results to a file
	print "\n\nOutputting sequences \n========================\n";
	open (OUTPUT, ">", "../output/processed/$family.txt") or die "Can't create output file\n";
	foreach my $keys (keys(%processed_results)){
		print OUTPUT ">$keys\n$processed_results{$keys}\n";
	}
	close (OUTPUT) or die "Can't close output file\n";
	print "\nData output to $family.txt\n";
	if (@ncbi_result_fail > 0) {
		open (FAILS, ">", "../output/$family-failed.txt") or die "Can't create output for failed searches\n";
		print FAILS "This file list sequences that were unable to be retrieved\n";
		foreach my $failed (@ncbi_result_fail) {
			print FAILS "$failed\n";
		}
		close FAILS or die "Can't close failed searches\n";
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
	