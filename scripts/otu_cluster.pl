#!/usr/bin/perl -w
# $Id: otu_cluster.pl 1 2008-12-03 22:42:29Z stajich $
use strict;

use Bio::AlignIO;

my %matrix;
my $needlefile = shift @ARGV;

my $in = Bio::AlignIO->new(-format => 'clustalw',
			   -file   => $needlefile);
while( my $aln = $in->next_aln ) {
    my ($q,$h) = $aln->each_seq;
    $matrix{$q->id}->{$q->id} = 0;
    next if $aln->no_sequences == 1;
    
    # this is only a pairwise alignment so I'm cheating here 
    # and pulling out only the first 2 seqs
    # average pairwise identity
    # $matrix{$q->id}->{$h->id} = sprintf("%.2f",$aln->average_percentage_identity);        
    my $pid = percent_id($h->seq, $q->seq);
    $matrix{$q->id}->{$h->id} = sprintf("%.8f", 1-$pid);
}

my @names = sort keys %matrix;
print scalar @names, "\n"; # print out the 1st line of DOTUR matrix which is
                           # number of names
for my $name ( @names ) {
    my @l;
    for my $n ( @names ) {
	push @l, $matrix{$name}->{$n};
    }
    print join("\t", $name, @l),"\n";
}

sub percent_id {
    my ($qseq,$hseq) = @_;
    if( length($qseq) ne length($hseq) ) {
	warn("cannot calculate PID if they aren't aligned!\n");
	return -1;
    }
    my $len = length($qseq);
    my $gaps = 0;
    my $match = 0;
    my $mismatch = 0;
    for(my $i = 0; $i < $len; $i++)  {
	my $qchar = lc substr($qseq,$i,1);
	my $hchar = lc substr($hseq,$i,1);
	if( $qchar eq '-' || $hchar eq '-' ||
	    $qchar eq '.' || $hchar eq '.' ) {
	    $gaps++;
	} elsif( $qchar eq $hchar ) {
	    $match++;
	} else {
	    $mismatch++;
	}
    }
    return $match / ($len - $gaps);
}
