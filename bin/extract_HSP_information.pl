############################################################
#  
# - a library of subroutines
#   from the examples and text in the book:
#
# Beginning Perl for Bioinformatics
# by James Tisdall
#
# published by O'Reilly & Associates
# (c) 2001 James Tisdall
#
# Version 20011230
#   incorporates a few errata and bug fixes
#
# Original code available in the "BeginPerlBioinfo.pm" module.
#
############################################################

#!/usr/bin/perl
use strict;
use warnings;

# Subroutine to extract HSP information
sub extract_HSP_information {
    my ($HSP) = @_;

    # declare and initialize variables
    my ($expect)      = '';
    my ($query)       = '';
    my ($query_range) = '';
    my ($subject)     = '';
    my ($subject_range) = '';

    ($expect) = ($HSP =~ /Expect = (\S+)/);

    $query = join( '', ($HSP =~ /^Query.*\n/gm) );

    $subject = join( '', ($HSP =~ /^Sbjct.*\n/gm) );

    my ($query_start, $query_end) = ($query =~ /(\d+).*\D(\d+)/s);
    my ($subject_start, $subject_end) = ($subject =~ /(\d+).*\D(\d+)/s);

    $query_range   = join( '..', $query_start, $query_end );
    $subject_range = join( '..', $subject_start, $subject_end );

    $query =~ s/[^acgt]//g;
    $subject =~ s/[^acgt]//g;

    return ( $expect, $query, $query_range, $subject, $subject_range, $query_start, $query_end );
}

# Main program
if ( @ARGV != 1 ) {
    die("Usage: $0 <blast_output_file>\n");
}

my $filename = $ARGV[0];

# Read the BLAST output file
open( my $fh, '<', $filename ) or die("Could not open file '$filename' $!\n");
my $blast_content = do { local $/; <$fh> };
close($fh);

# Split the content into alignments based on the sequence identifier
my @alignments = split(/^BLASTN.*?^Query=.*?^>.*?\n/ms, $blast_content);

# Process each alignment
my @HSP_positions; # to collect start and end positions

foreach my $alignment (@alignments) {
    next unless $alignment =~ /\S/;  # Skip empty alignments

    # Parse HSPs
    my @HSPs = parse_blast_alignment_HSP($alignment);

    # Process each HSP
    foreach my $hsp (@HSPs) {
        my ( $expect, $query, $query_range, $subject, $subject_range, $query_start, $query_end ) =
          extract_HSP_information($hsp);

        # Print the results for each HSP
        print "\n-> Expect value:   $expect\n";
        print "\n-> Query string:   $query\n";
        print "\n-> Query range:    $query_range\n";
        print "\n-> Subject String: $subject\n";
        print "\n-> Subject range:  $subject_range\n";

        # Print the start and end positions
        print "\n-> Query start position: $query_start\n";
        print "\n-> Query end position:   $query_end\n";

        # Collect the start and end positions
        push @HSP_positions, [$query_start, $query_end];
    }
}

# Determine the overall range on the query sequence
my @sorted_positions = sort { $a <=> $b } map { $_->[0], $_->[1] } @HSP_positions;
my $overall_start = $sorted_positions[0];
my $overall_end   = $sorted_positions[-1];

print "\nOverall range on query sequence: $overall_start..$overall_end\n";

