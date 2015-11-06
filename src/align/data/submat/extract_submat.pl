#!/usr/bin/env perl

# extract a hard-coded substitution matrix from a source file of BLAST+
# e.g. ncbi-blast-2.2.31+-src/c++/src/util/tables/sm_blosum62.c

use strict;
use warnings;
use File::Basename;

print "# filename: ", (basename $ARGV[0]), "\n";
while (<>) {
    if (/static const TNCBIScore (\w+)/) {
        print "# variable: $1\n";
        my $header = <> . <>;
        my @aa = ($header =~ /([A-Z*])/g);
        # remove the first and last *
        shift @aa; pop @aa;
        print " ";
        for (@aa) {
            printf "%4s", $_;
        }
        print "\n";
    }
    elsif (/\/\*([A-Z*])\*\//) {
        my $scores = $_ . <>;
        print $1;
        my @scores = ($scores =~ /(-?\d+)/g);
        for (@scores) {
            printf "%4d", $_;
        }
        print "\n";
    }
}
