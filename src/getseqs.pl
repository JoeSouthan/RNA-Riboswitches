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
use lib 'Modules';
use GetSeqs; 
use ProcessSeqs;
use Notices;
#
#	Script Logic
#

	my ($mode, @rfam_families);
	system("cls");
	#Turn off buffer for status messages
		$|=1;
	#Check for params
		Error(0) unless (@ARGV);
	#Default to seed
		$mode = ($ARGV[0] ne "seed" or $ARGV[0] ne "full") ? $ARGV[0] : "seed" ;
	#Ensure search file exists
		Error(1) unless (defined($ARGV[1]));
		@rfam_families = @{OpenSearch($ARGV[1])};
	#Warn the user
		if ($mode eq "full"){
			Note(0);
			my $choice = <STDIN>;
			chomp ($choice);
			$mode = ($choice eq "n") ? "seed" : "full";
		}

#
#	Main loop
#

	my $rfam_families_count = scalar(@rfam_families);
	foreach my $family (@rfam_families) {
		system("cls");
		print "\n========================\nFamily: $family $rfam_families_count remain(s).\n========================\n";

		#
		#	Fetch RFAM
		#
			print "\n\nRetrieving $mode alignments from RFAM\n========================\n";
				my @rfam_result_ref = FetchRFAM($family, $mode);
				my %rfam_result = %{$rfam_result_ref[0]};
				my @gis = keys (%rfam_result);
	 			my $no_ids = scalar(@gis);
			print "Found $rfam_result_ref[1] records from RFAM for $family family. ($no_ids individual genes)\n";

		#
		#	Fetch NCBI
		#
			print "\n\nRetrieving sequences from NCBI for $family \n========================\n\r";
				my $ncbi_result = FetchNCBI($family, \@gis);

		#
		#	Split GenBank file
		#
			if (0 == $ncbi_result) {
				Note(1);
				SplitGenBank($family, $no_ids);
			}

		#
		# Process the sequences
		#
			ProcessSequences($family, \%rfam_result);

		$rfam_families_count--;
		sleep 3;
	}
