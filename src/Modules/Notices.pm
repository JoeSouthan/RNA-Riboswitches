#! /usr/bin/perl -w
#
# Notices.pm
# Created by: Joseph Southan
# Date: 27/6/13
# Description: - Module for getseqs.pl
# 
package Notices;
use strict;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw (Error Note);

sub Error {
	my @errors = (
		#Logic
		"\n========================\nError: Please ensure to have the correct arguments\n========================\n",
		"\n========================\nError: No search file found\n========================\n",
		#GetSeq
		"Can't open search file, ensure it exists.\n",
		"Can't close search file.\n"
	);
	die $errors[$_[0]];
}
sub Note {
	my @notes = (
		#getseqs.pl
		"\n========================\nAre you sure you want to get the full alignments? \nThis is very data and processor hungry and some sequences may not contain the full alignments. y/n \n",
		"\n\nSplitting GenBank file \n========================\n\r",
		#GetSeqs.pm
		"Getting sequences, please wait...\r",
		"Raw files for this family exist. Overwrite?\t"
	);
	print $notes[$_[0]];
}
1;
