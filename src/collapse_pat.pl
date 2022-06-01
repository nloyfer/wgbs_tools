#!/usr/bin/perl

=for

This scripts takes a sorted (-k2,2n -k3,3) uncompressed pat file,
with possibly "almost duplicated" lines, and collapses them as follows:
For example, there are adjacent lines in the input file:
chr1    90  CCC 2
chr1    90  CCC 3
Then it will output:
chr1    90  CCC 5

The resulted pat is printed to stdout

=cut

use strict;
use warnings;

my $file = $ARGV[0];
if (not defined $file) {
    die "Need file name\n";
}

# set pre_line
my @pre_line;

# read whole file
my $count = 0;
open my $info, $file or die "Could not open $file: $!";

while( my $line = <$info>)  {
    chomp $line;
    my @words = split(/\t/, $line);

    # First line init
    if ($. == 1) {
	@pre_line = @words;
	$pre_line[3] = 0;
    }

    # compare current line to previous one
    my $equals = 1;

    foreach (my $i = 0; $i < @words; $i++) {
	
        if ($words[$i] ne $pre_line[$i] && $i != 3) {
            $equals = 0;
            last;
        }
    }
    if ($equals) {
        $count += $words[3];
    }
    else {
	if ($count > 0) {
 	        $pre_line[3] = $count;
        	print join( "\t", @pre_line ), "\n";
	}
	$count = $words[3];
    }
    @pre_line = @words
}

# print last line:
if ($count > 0) {
	$pre_line[3] = $count;
	print join( "\t", @pre_line ), "\n";
}

close $info;
