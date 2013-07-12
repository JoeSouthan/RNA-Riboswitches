#! /usr/bin/perl -w
#
# GetSeqs.pm
# Created by: Joseph Southan
# Date: 27/6/13
# Description: - Module for getseqs.pl
# 
package GetSeqs;
use strict;
use LWP::Simple;
use Data::Dumper;
use Bio::SeqIO;
use LWP::UserAgent;
use File::Path qw(make_path remove_tree);
use Notices;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw (FetchRFAM FetchNCBI OpenSearch RFAM_ErrorOut);

sub OpenSearch {
	my ($search_file) = @_;
	my @rfam_families;
	#Open up the search file
	open (SEARCHES, "<", $search_file) or Error(2);
	while (my $line = <SEARCHES>){
		chomp ($line);
		push (@rfam_families, $line) if (length ($line) > 2);
	}
	close (SEARCHES) or Error(3);
	return \@rfam_families;
}

#$family
#returns hash ref, $count
sub FetchRFAM {
	#Read in the RFAM list
	#http://rfam.sanger.ac.uk/family/RF00174/alignment/seed/stockholm?alnType=seed&nseLabels=0&cs=0
	my ($family, $mode, $search_file) = @_;
	my (%rfam_result, @rfam_families);
	my ($count, $exist_count) = (0, 0);

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
				push( @{$rfam_result{$1}}, @sub_result); 
				$exist_count++;
			} else {
				my @sub_result = [@locs];
				$rfam_result{$1} = [@sub_result];
			}
			$count++;
		}
		return (\%rfam_result, $count);
	} else {
		print "$family Get Failed\n";
		$rfam_result{"error"} = 1;
		return (\%rfam_result, 0);
	}
}

sub FetchNCBI {
	my $family = $_[0];
	my @gis = @{$_[1]};
	my $existing_raw = "../output/raw/$family.txt";
	my $existing_raw_split ="../output/raw/$family";
	if (-e $existing_raw or -d $existing_raw_split) {
		Note(3);
		my $overwrite_choice = <STDIN>;
		chomp $overwrite_choice;
		if ($overwrite_choice eq "y") {
			unlink($existing_raw);
			if (-d $existing_raw_split) {
				remove_tree($existing_raw_split);
			}
			return _FetchSeqs(\@gis, $family);
		} else {
			return 2; #Use existing files
		}
	} else {
		return _FetchSeqs(\@gis, $family);
	}
}
sub _FetchSeqs {
	my @gis = @{$_[0]};
 	my $family = $_[1];
 	my $ua = LWP::UserAgent->new();
 	my $fetch_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gbwithparts&id=";

 	$fetch_url .= join (",", @gis);

 	#Progres bar
 	#$ua->show_progress('true value');
 	Note(2);
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
sub RFAM_ErrorOut {
	open (ERROR, ">>", "RFAM_Failed.txt") or die ("Can't open error output\n");
	print ERROR localtime().": $_[0] failed to download\n";
	close ERROR or die ("Can't close error output\n");
}
1;
