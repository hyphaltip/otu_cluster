#!/usr/bin/perl -w
# $Id: needle2dotur.pl 1 2008-12-03 22:42:29Z stajich $
# Usage: needle2dotur.pl report.needle

# Author: Jason Stajich <jason[_AT_]open-bio.org>
# http://fungalgenomes.org/

# this requires the all-vs-all needleman-wunsch report should have been already run
# see run_do_needle.sh for script to run all the pairwise alignments using
# needle (part of the EMBOSS package) http://www.emboss.org/

use strict;
use Getopt::Long;


my $Usage = <<EOF
needle_clusters.pl [options] <emboss report> > report.doturmatrix 
    Options
    -h     : show brief help
    -o     : dotur distance matrix filename (defaults to STDOUT)
    -d     : distance instead of percent identity
EOF
    ;							

my ($Distance,$outfile);
GetOptions(
	   'h|help'      => sub { print STDERR $Usage,"\n";
				  exit(0) },
	   'o|out:s'     => \$outfile,
	   'd|distance!' => \$Distance,
	   );
my $ofh;
if( $outfile ) {
    open($ofh, ">$outfile") || die "cannot open $outfile for writing: $!";
} else {
    $ofh = \*STDOUT; # default to STDOUT for writing report
}

my ($reportfile) = shift @ARGV;  # get the cmdline remaining arguments
				    # from ARGV

if( ! defined $reportfile ) { 
    print STDERR $Usage, "\n";
    exit(0);
}


my %matrix    = &read_emboss_identity($reportfile);

my @names = sort keys %matrix;
print scalar @names, "\n"; # print out the 1st line of DOTUR matrix which is
                           # number of names
for my $name ( @names ) {
    my @l;
    for my $n ( @names ) {
	if( ! exists $matrix{$name}->{$n} ) {
	    warn("cannot find $name->$n\n");
	    next;
	}
	push @l, $matrix{$name}->{$n};
    }
    print join("\t", $name, @l),"\n";
}


# --- Functions below ---

# parse Emboss needle or water output
# return a hash with key being the sequence id,
#         each value is in tern a hashref pointing to the sequences
#         which were aligned and the value for this $hash{key}->{key} will be
#         the percent identity 
# input is the name of the file
sub read_emboss_identity {
    my ($emboss_report,$no_skip_same) = @_;
    open(REPORT, "cat $emboss_report |") || die("cannot open file $emboss_report: $!");
    my (@seqs,%lengths,%gaps,%identities,%distmat);
    while(<REPORT>) {
	if( /^\#\s+Aligned_sequences/) {
	    @seqs = ();
	} elsif( /^\#\s+(1|2):\s+(\S+)/) {
	    $seqs[$1-1] = $2;
	} elsif( /^\#\s+Length:\s+(\d+)/ ) {
	    $lengths{$seqs[0]}->{$seqs[1]} = $1;
	} elsif( /^\#\s+Identity:\s+(\d+)\/(\d+)/ ) {
	    $identities{$seqs[0]}->{$seqs[1]} = $1;	    
	} elsif( /^\#\s+Gaps:\s+(\d+)\/(\d+)/ ) {
	    $gaps{$seqs[0]}->{$seqs[1]} = $1;
	}
    }
    close(REPORT);
    for my $seqa ( keys %lengths ) {
	for  my $seqb ( keys %{$lengths{$seqa}} ) {
	    my $id = ( $identities{$seqa}->{$seqb} / 
		       ($lengths{$seqa}->{$seqb} - $gaps{$seqa}->{$seqb}));
	    if( $Distance  ) {
		$id = 1 - $id;
	    }
	    $distmat{$seqa}->{$seqb} = sprintf("%.8f",$id);
	}
    }
    return %distmat;
}
