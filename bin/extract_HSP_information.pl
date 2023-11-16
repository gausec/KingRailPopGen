############################################################
#
# BeginPerlBioinfo.pm
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

    $query_range   = join( '..', ($query =~ /(\d+).*\D(\d+)/s) );
    $subject_range = join( '..', ($subject =~ /(\d+).*\D(\d+)/s) );

    $query =~ s/[^acgt]//g;
    $subject =~ s/[^acgt]//g;

    return ( $expect, $query, $query_range, $subject, $subject_range );
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
foreach my $alignment (@alignments) {
    next unless $alignment =~ /\S/;  # Skip empty alignments

    # Parse HSPs
    my @HSPs = parse_blast_alignment_HSP($alignment);

    # Process each HSP
    foreach my $hsp (@HSPs) {
        my ( $expect, $query, $query_range, $subject, $subject_range ) =
          extract_HSP_information($hsp);

        # Print the results for each HSP
        print "\n-> Expect value:   $expect\n";
        print "\n-> Query string:   $query\n";
        print "\n-> Query range:    $query_range\n";
        print "\n-> Subject String: $subject\n";
        print "\n-> Subject range:  $subject_range\n";
    }
}

# Subroutine to parse HSPs from a BLAST alignment
sub parse_blast_alignment_HSP {
    my ($alignment) = @_;

    # declare and initialize variables
    my $beginning_annotation = '';
    my $HSP_section          = '';
    my @HSPs                 = ();

    # Extract the beginning annotation and HSPs
    ( $beginning_annotation, $HSP_section ) =
      ( $alignment =~ /(.*?)(^ Score =.*)/ms );

    # Store the $beginning_annotation as the first entry in @HSPs
    push( @HSPs, $beginning_annotation );

    # Parse the HSPs, store each HSP as an element in @HSPs
    while ( $HSP_section =~ /(^ Score =.*\n)(^(?! Score =).*\n)+/gm ) {
        push( @HSPs, $& );
    }

    # Return an array with first element = the beginning annotation,
    # and each successive element = an HSP
    return @HSPs;
}
